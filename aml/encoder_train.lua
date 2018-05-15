-- Train an encoder using config.json.

require('torch')
require('nn')
require('optim')
require('hdf5')
require('cunn')
require('lfs')

package.path = package.path .. ";" .. lfs.currentdir() .. '/../?/th/init.lua'
lib = require('lib')

--- Check if the file exists and read it.
-- @param filepath path to H5 file
-- @return tensor read from file
local function checkAndRead(filepath)
  if not lib.utils.fileExists(filepath) then
    print('[Error] file ' .. filepath .. ' not found; check ' .. configFile)
    os.exit()
  end

  print('[Training] reading ' .. filepath)
  return lib.utils.readHDF5(filepath)
end

-- https://github.com/torch/cutorch
if cutorch.isCachingAllocatorEnabled() then
  print('[Training] caching allocator enabled')
end

--- Convert a table of tensors to a table of CUDA tensors.
-- @param table the table to convert
local function tableToCUDA(table)
  for key, value in pairs(table) do
    table[key] = table[key]:cuda()
  end
end

--- Set all table elements to zero.
-- @param table table to set to zero
local function tableToZero(table)
  for key, value in pairs(table) do
    table[key] = 0
  end
end

-- ---------------------------------------------------------------
-- Read configuration and data.
-- ---------------------------------------------------------------

-- Load configuration.
configFile = 'clean.json'
if arg[1] then
  configFile = arg[1]
end

if not lib.utils.fileExists(configFile) then
  print('[Error] configuration file not found')
  os.exit()
end

print('[Training] reading ' .. configFile)
config = lib.utils.readJSON(configFile)

config['data_directory'] = config['data_directory'] .. '/'
config['base_directory'] = config['base_directory'] .. '/'

if not lib.utils.directoryExists(config['data_directory']) then
  print('[Error] Data directory not found')
  os.exit()
end

if not lib.utils.fileExists(config['base_directory'] .. config['prior_model_file']) then
  print('[Error] ' .. config['base_directory'] ..  config['prior_model_file'] .. ' does not exist; train a shape prior first, create the base directory manually and put the trained model inside')
  os.exit()
end

if lib.utils.fileExists(config['base_directory'] .. config['inference_model_file']) then
  print('[Training] Model already exists, overwriting results and model?')
  print('[Training] Press ANY KEY to continue')
  io.read()
end

config['data_directory'] = config['data_directory'] .. '/'
points = checkAndRead(config['data_directory'] .. config['inference_training_inputs'])
inputs_2 = checkAndRead(config['data_directory'] .. config['inference_training_lsdf_inputs'])
inputs = torch.cat({points:float(), inputs_2}, 2)

if config['inference_training_outputs'] and config['inference_training_lsdf_outputs'] then
  outputs_1 = checkAndRead(config['data_directory'] .. config['inference_training_outputs'])
  outputs_2 = checkAndRead(config['data_directory'] .. config['inference_training_lsdf_outputs'])
  outputs = torch.cat({outputs_1:float(), outputs_2}, 2)
end

space = checkAndRead(config['data_directory'] .. config['inference_training_space'])

-- Load data for validation.
valPoints = checkAndRead(config['data_directory'] .. config['validation_inputs'])
valInputs_2 = checkAndRead(config['data_directory'] .. config['validation_lsdf_inputs'])
valInputs = torch.cat({valPoints:float(), valInputs_2}, 2)

if config['validation_outputs'] and config['validation_lsdf_outputs'] then
  valOutputs_1 = checkAndRead(config['data_directory'] .. config['validation_outputs'])
  valOutputs_2 = checkAndRead(config['data_directory'] .. config['validation_lsdf_outputs'])
  valOutputs = torch.cat({valOutputs_1:float(), valOutputs_2}, 2)
end

valSpace = checkAndRead(config['data_directory'] .. config['validation_space'])

-- Check dimensions.
N = inputs:size(1)

-- ---------------------------------------------------------------
-- Setup and initialize model.
-- ---------------------------------------------------------------

priorModel = torch.load(config['base_directory'] .. config['prior_model_file'])
print(priorModel)

encoder = priorModel.modules[1]
encoder:remove(#encoder.modules)
encoder:remove(#encoder.modules)
mean = encoder.modules[#encoder.modules].modules[1]
encoder:remove(#encoder.modules)
encoder:add(mean)
decoder = priorModel.modules[2]

KLD = nn.MaximumLikelihoodPrior()
KLD.var = 1
KLD.lambda = config['prior_weight']
encoder:add(KLD)

lib.utils.fixLayersAfter(decoder, 1)
if config['reinitialize_encoder'] then
  lib.init(encoder, config['weight_initialization'], config['weight_value'], config['bias_initialization'], config['bias_value'])
  print('[Training] reinitializing encoder')
end

model = nn.Sequential()
model:add(encoder)
model:add(decoder)

model = model:cuda()
print(model)

-- ---------------------------------------------------------------
-- Setup training.
-- ---------------------------------------------------------------

batchSize = config['batch_size']
learningRate = config['learning_rate']
momentum = config['momentum']
weightDecay = config['weight_decay']
epochs = config['epochs']
iterations = epochs*math.floor(N/batchSize)
minimumLearningRate = config['minimum_learning_rate']
learningRateDecay = config['decay_learning_rate']
maximumMomentum = config['maximum_momentum']
momentumDecay = config['decay_momentum']
decayIterations = config['decay_iterations']
lossIterations = config['loss_iterations']
snapshotIterations = config['snapshot_iterations']
testIterations = config['test_iterations']

if config['inference_training_statistics'] then
  if not lib.utils.fileExists(config['data_directory'] .. config['inference_training_statistics']) then
    print('[Error] the key inference_training_statistics is defined in the configuration; but file ' .. config['data_directory'] ..config['inference_training_statistics'] .. ' cannot be found')
    os.exit()
  end

  print('[Training] reading ' .. config['data_directory'] ..config['inference_training_statistics'])
  trainStatistics = lib.utils.readHDF5(config['data_directory'] ..config['inference_training_statistics'])
  trainStatistics = nn.utils.addSingletonDimension(nn.utils.addSingletonDimension(trainStatistics, 1), 1)
  trainStatistics = torch.repeatTensor(trainStatistics, batchSize, 1, 1, 1, 1)
end

-- Parameters and gradients.
parameters, gradParameters = model:getParameters()
parameters = parameters:cuda()
gradParameters = gradParameters:cuda()

-- ---------------------------------------------------------------
-- Setup loss.
-- ---------------------------------------------------------------

criterion = nn.MultipleCriterion()
criterion.weights = {}
criterion.criteria = {}
criterion.channels = {}

assert(#config['criteria'] == #config['weights'])
for i = 1, #config['criteria'] do
  if config['criteria'][i] == 'sdfpointbce' then
    local pointCriterion = nn.GaussianPointBCECriterion()
    pointCriterion.sizeAverage = false
    pointCriterion = pointCriterion:cuda()

    criterion.criteria[i] = pointCriterion
    criterion.weights[i] = config['weights'][i]
    criterion.channels[i] = 2

  elseif config['criteria'][i] == 'sdfspacebce' then
    local spaceCriterion = nn.GaussianFreeSpaceBCECriterion()
    spaceCriterion.sizeAverage = false
    spaceCriterion = spaceCriterion:cuda()

    if config['weighted'] then
      assert(trainStatistics)

      spaceCriterion = nn.WeightedGaussianFreeSpaceBCECriterion()
      spaceCriterion.sizeAverage = false
      spaceCriterion.weights = trainStatistics
      spaceCriterion = spaceCriterion:cuda()

      print('[Training] using weighted sdfspacebce')
    end

    criterion.criteria[i] = spaceCriterion
    criterion.weights[i] = config['weights'][i]
    criterion.channels[i] = 2

  elseif config['criteria'][i] == 'occpointbce' then
    local pointCriterion = nn.PointBCECriterion()
    pointCriterion.sizeAverage = false
    pointCriterion = pointCriterion:cuda()

    criterion.criteria[i] = pointCriterion
    criterion.weights[i] = config['weights'][i]
    criterion.channels[i] = 1

  elseif config['criteria'][i] == 'occspacebce' then
    local spaceCriterion = nn.FreeSpaceBCECriterion()
    spaceCriterion.sizeAverage = false
    spaceCriterion = spaceCriterion:cuda()

    if config['weighted'] then
      assert(trainStatistics)

      spaceCriterion = nn.WeightedFreeSpaceBCECriterion()
      spaceCriterion.sizeAverage = false
      spaceCriterion.weights = trainStatistics
      spaceCriterion = spaceCriterion:cuda()

      print('[Training] using weighted occspacebce')
    end

    criterion.criteria[i] = spaceCriterion
    criterion.weights[i] = config['weights'][i]
    criterion.channels[i] = 1

  else
    assert(False)
  end
end

-- ---------------------------------------------------------------
-- Main training loop.
-- ---------------------------------------------------------------

protocol = torch.Tensor(iterations, 5):fill(0)
valProtocol = torch.Tensor(math.floor(iterations/testIterations) + 1, 10):fill(0)

for t = 1, iterations do

  -- Sample a random batch from the dataset.
  local shuffle = torch.randperm(N)
  shuffle = shuffle:narrow(1, 1, batchSize)
  shuffle = shuffle:long()

  local container = {} -- Container will hold the different inputs and outputs in different representations.
  container['input'] = inputs:index(1, shuffle)
  container['points'] = points:index(1, shuffle)
  container['space'] = space:index(1, shuffle)

  if outputs then
    container['target'] = outputs:index(1, shuffle)
  end

  tableToCUDA(container)

  --- Definition of the objective on the current mini-batch.
  -- This will be the objective fed to the optimization algorithm.
  -- @param x input parameters
  -- @return object value, gradients
  local feval = function(x)

    if x ~= parameters then
      parameters:copy(x)
    end

    gradParameters:zero()

    -- Evaluate function on mini-batch.
    local pred = model:forward(container['input'])
    local f = criterion:forward(pred, {container['points'], container['space'], container['points'], container['space']})
    local df = criterion:backward(pred, {container['points'], container['space'], container['points'], container['space']})
    model:backward(container['input'], df)

    -- Save losses in protocol.
    protocol[t][2] = f
    protocol[t][3] = KLD.loss

    -- Weight decay:
    if weightDecay > 0 then
      weightDecayLoss = weightDecay * torch.norm(parameters, 2)^2/2
      f = f + weightDecayLoss

      protocol[t][4] = weightDecayLoss
      gradParameters:add(parameters:clone():mul(weightDecay))
    end

    if container['target'] then
      protocol[t][5] = torch.mean(torch.abs(pred - container['target']))
    end

    -- Add the Kullback-Leibler loss.
    f = f + KLD.loss

    return f, gradParameters
  end

  -- Update state with learning rate and momentum.
  state = state or {
    learningRate = learningRate,
    momentum = momentum,
    learningRateDecay = 0 -- will be done manually below
  }

  -- Returns the new parameters and the objective evaluated
  -- before the update.
  --_, _ = optim.adam(feval, parameters, state)
  _, _ = optim.sgd(feval, parameters, state)
  local time = os.date("*t")

  -- Report a smoothed loss instead of batch loss.
  if t%lossIterations == 0 then
    local smoothedLoss = torch.mean(protocol:narrow(1, t - lossIterations + 1, lossIterations):narrow(2, 2, 1))
    local smoothedKLD = torch.mean(protocol:narrow(1, t - lossIterations + 1, lossIterations):narrow(2, 3, 1))
    local smoothedWeightDecay = torch.mean(protocol:narrow(1, t - lossIterations + 1, lossIterations):narrow(2, 4, 1))
    local loss = smoothedLoss + smoothedKLD + smoothedWeightDecay

    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec ..': [' .. KLD.lambda .. '] ' .. loss .. ' / ' .. smoothedLoss .. ' / ' .. smoothedKLD .. ' / ' .. smoothedWeightDecay)
  end

  -- Decay learning rate and KLD weight, do this before resetting all smoothed
  -- statistics.
  if t%decayIterations == 0 then
    learningRate = math.max(minimumLearningRate, learningRate*learningRateDecay)
    momentum = math.min(maximumMomentum, momentum*momentumDecay)
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec ..': learning rate ' .. learningRate)
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec ..': momentum ' .. momentum)
  end

  -- Validate on validation set.
  if t%testIterations == 0 or t == 1 or t == iterations then

    -- In case the validation set gets to large.
    local valN = valInputs:size(1)
    local valBatchSize = batchSize
    local valNumBatches = math.floor(valN/valBatchSize)
    local valIteration = math.floor(t/testIterations) + 1

    -- Accumulate and save all predictions.
    local accValPreds = nil

    for b = 0, valNumBatches do
      local batchStart = b*valBatchSize + 1
      local batchLength = math.min((b + 1)*valBatchSize - b*valBatchSize, valN - b*valBatchSize)

      local container = {} -- Container will hold the different inputs and outputs in different representations.
      container['input'] = valInputs:narrow(1, batchStart, batchLength)
      container['points'] = valPoints:narrow(1, batchStart, batchLength)
      container['space'] = valSpace:narrow(1, batchStart, batchLength)

      if valTargets then
        container['target'] = valTargets:narrow(1, batchStart, batchLength)
      end

      tableToCUDA(container)

      local valPreds = model:forward(container['input'])
      local valCodes = model.modules[1].modules[#model.modules[1].modules].output
      local f = criterion:forward(valPreds, {container['points'], container['space'], container['points'], container['space']})

      accValPreds = lib.utils.appendTensor(accValPreds, valPreds)

      valProtocol[valIteration][2] = valProtocol[valIteration][2] + f
      valProtocol[valIteration][4] = valProtocol[valIteration][4] + KLD.loss

      if container['target'] then
        valProtocol[valIteration][6] = valProtocol[valIteration][6] + torch.mean(torch.abs(valPreds - container['target']))
      end

      valProtocol[valIteration][8] = valProtocol[valIteration][8] + torch.mean(valCodes)
      valProtocol[valIteration][9] = valProtocol[valIteration][9] + torch.std(valCodes)
    end

    valProtocol[valIteration][1] = t
    for i = 2, valProtocol:size(2) do
      valProtocol[valIteration][i] = valProtocol[valIteration][i] / valNumBatches
    end

    local predFile = config['base_directory'] .. t .. '_predictions.h5'
    lib.utils.writeHDF5(predFile, accValPreds)
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec ..': saved ' .. predFile)

    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec ..': codes ' .. valProtocol[valIteration][8] .. ' / ' .. valProtocol[valIteration][9])
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec ..': validation loss ' .. valProtocol[valIteration][2])
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec ..': KLD loss ' .. valProtocol[valIteration][4])
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec ..': validation error (abs) ' .. valProtocol[valIteration][6])
  end
end

-- ---------------------------------------------------------------
-- Save protocols and model.
-- ---------------------------------------------------------------

protocolFile = config['base_directory'] .. config['train_protocol_file']
lib.utils.writeHDF5(protocolFile, protocol)
print('[Training] protocol ' .. protocolFile)

protocolFile = config['base_directory'] .. config['val_protocol_file']
lib.utils.writeHDF5(protocolFile, valProtocol)
print('[Training] protocol ' .. protocolFile)

modelFile = config['base_directory'] .. config['inference_model_file']
torch.save(modelFile, model)
print('[Training] model ' .. modelFile)