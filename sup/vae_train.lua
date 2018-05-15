-- Train an auto-encoder using config.json.

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

if lib.utils.directoryExists(config['base_directory']) then
  print('[Training] Base directory already exists, overwriting results and model?')
  print('[Training] Press ANY KEY to continue')
  io.read()
end

inputs_1 = checkAndRead(config['data_directory'] .. config['inference_training_inputs'])
inputs_2 = checkAndRead(config['data_directory'] .. config['inference_training_lsdf_inputs'])
inputs = torch.cat({inputs_1:float(), inputs_2}, 2)

outputs_1 = checkAndRead(config['data_directory'] .. config['inference_training_outputs'])
outputs_2 = checkAndRead(config['data_directory'] .. config['inference_training_lsdf_outputs'])
outputs = torch.cat({outputs_1:float(), outputs_2}, 2)

-- Load data for validation.
valInputs_1 = checkAndRead(config['data_directory'] .. config['validation_inputs'])
valInputs_2 = checkAndRead(config['data_directory'] .. config['validation_lsdf_inputs'])
valInputs = torch.cat({valInputs_1:float(), valInputs_2}, 2)

valOutputs_1 = checkAndRead(config['data_directory'] .. config['validation_outputs'])
valOutputs_2 = checkAndRead(config['data_directory'] .. config['validation_lsdf_outputs'])
valOutputs = torch.cat({valOutputs_1:float(), valOutputs_2}, 2)

-- Create snapshot directory.
print('[Training] creating ' .. config['base_directory'])
if not lib.utils.directoryExists(config['base_directory']) then
  lib.utils.makeDirectory(config['base_directory'])
end

-- ---------------------------------------------------------------
-- Setup and initialize model.
-- ---------------------------------------------------------------

-- For later simplicity.
N = inputs:size()[1]
channels = inputs:size()[2]
height = inputs:size()[3]
width = inputs:size()[4]
depth = inputs:size()[5]

-- Set up config for model.
autoEncoderConfig = lib.variationalAutoEncoder.config
autoEncoderConfig.height = height
autoEncoderConfig.width = width
autoEncoderConfig.depth = depth
autoEncoderConfig.code = config['code']
autoEncoderConfig.channels = config['channels']
autoEncoderConfig.kernelSizes = config['kernel_sizes']
autoEncoderConfig.pooling = config['pooling']
autoEncoderConfig.poolingSizes = config['pooling_sizes']
autoEncoderConfig.transfer = nn.ReLU
autoEncoderConfig.transfers = config['transfers']
autoEncoderConfig.normalizations = config['normalizations']
autoEncoderConfig.dataChannels = channels
autoEncoderConfig.outputChannels = channels
--autoEncoderConfig.printDimensions = true
--autoEncoderConfig.checkNaN = true

-- Set up the auto-encoder.
model = nn.Sequential()
model, context = lib.variationalAutoEncoder.autoEncoder(model, autoEncoderConfig)

nonLinearities = nn.PerChannelNonLinearity()
nonLinearities.layers = {nn.Sigmoid(), nn.Identity()}
context['decoder']:add(nonLinearities)

KLD = context['KLD']
KLD.sizeAverage = false
KLD.lambda = config['prior_weight']
print('[Training] using prior weight ' .. KLD.lambda)

mean = context['mean']
logvar = context['logVar']
decoder = context['decoder']

-- Initialize weights.
print('[Training] initializing model')
lib.init(model, config['weight_initialization'], config['weight_value'],
  config['bias_initialization'], config['bias_value'])
model = model:cuda()
print(model)

-- ---------------------------------------------------------------
-- Setup loss.
-- ---------------------------------------------------------------

-- Criterion.
occCriterion = nn.BCECriterion()
occCriterion.sizeAverage = false
occCriterion = occCriterion:cuda()

sdfCriterion = nn.FixedVarianceGaussianNLLCriterion()
sdfCriterion.sizeAverage = false
sdfCriterion.logvar = 0
sdfCriterion = sdfCriterion:cuda()

criterion = nn.PerChannelCriterion()
criterion.criteria = {occCriterion, sdfCriterion }
criterion.weights = {1, 1 }
criterion = criterion:cuda()

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
testIterations = config['test_iterations']

parameters, gradParameters = model:getParameters()
parameters = parameters:cuda()
gradParameters = gradParameters:cuda()

-- ---------------------------------------------------------------
-- Main training loop.
-- ---------------------------------------------------------------

protocol = torch.Tensor(iterations, 5):fill(0)
valProtocol = torch.Tensor(math.floor(iterations/testIterations) + 1, 10):fill(0)

-- Main training loop.
for t = 1, iterations do

  -- Sample a random batch from the dataset.
  local shuffle = torch.randperm(N)
  shuffle = shuffle:narrow(1, 1, batchSize)
  shuffle = shuffle:long()

  local input = inputs:index(1, shuffle)
  local output = outputs:index(1, shuffle)

  input = input:cuda()
  output = output:cuda()
  protocol[t][1] = t

  --- Definition of the objective on the current mini-batch.
  -- This will be the objective fed to the optimization algorithm.
  -- @param x input parameters
  -- @return object value, gradients
  local feval = function(x)

    -- Get new parameters.
    if x ~= parameters then
      parameters:copy(x)
    end

    -- Reset gradients
    gradParameters:zero()

    -- Evaluate function on mini-batch.
    local pred = model:forward(input)
    local f = criterion:forward(pred, output)

    protocol[t][2] = f
    protocol[t][3] = KLD.loss

    -- Estimate df/dW.
    local df_do = criterion:backward(pred, output)
    model:backward(input, df_do)

    -- Weight decay:
    if weightDecay > 0 then
      weightDecayLoss = weightDecay * torch.norm(parameters,2)^2/2
      protocol[t][4] = weightDecayLoss

      f = f + weightDecayLoss
      gradParameters:add(parameters:clone():mul(weightDecay))
    end

    protocol[t][5] = torch.mean(torch.abs(pred - output))

    -- Add the Kullback-Leibler loss.
    f = f + KLD.loss

    -- return f and df/dX
    return f, gradParameters
  end

  -- Update state with learning rate and momentum.
  adamState = adamState or {
    learningRate = learningRate,
    momentum = momentum,
    learningRateDecay = 0 -- will be done manually below
  }

  -- Returns the new parameters and the objective evaluated
  -- before the update.
  --p, f = optim.adam(feval, parameters, adamState)
  p, f = optim.sgd(feval, parameters, adamState)
  local time = os.date("*t")

  -- Report a smoothed loss instead of batch loss.
  if t%lossIterations == 0 then
    local smoothedLoss = torch.mean(protocol:narrow(1, t - lossIterations + 1, lossIterations):narrow(2, 2, 1))
    local smoothedKLD = torch.mean(protocol:narrow(1, t - lossIterations + 1, lossIterations):narrow(2, 3, 1))
    local smoothedWeightDecay = torch.mean(protocol:narrow(1, t - lossIterations + 1, lossIterations):narrow(2, 4, 1))
    local loss = smoothedLoss + smoothedKLD + smoothedWeightDecay

    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': ' .. loss .. ' / ' .. smoothedLoss .. ' / ' .. smoothedKLD .. ' / ' .. smoothedWeightDecay)
  end

  -- Decay learning rate and KLD weight, do this before resetting all smoothed
  -- statistics.
  if t%decayIterations == 0 then
    learningRate = math.max(minimumLearningRate, learningRate*learningRateDecay)
    momentum = math.min(maximumMomentum, momentum*momentumDecay)
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': learning rate ' .. learningRate)
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': momentum ' .. momentum)
  end

  -- Validate on validation set.
  if t%testIterations == 0 or t == 1 or t == iterations then

    -- In case the validation set gets to large.
    local valN = valInputs:size(1)
    local valBatchSize = batchSize
    local valNumBatches = math.floor(valN/valBatchSize)
    local valIteration = math.floor(t/testIterations) + 1

    local accValMeanPreds = nil

    -- iteration, loss x2, KLD loss x2, error x2, thresh error x2, mean mean, var mean, mean logvar
    for b = 0, valNumBatches do
      local batchStart = b*valBatchSize + 1
      local batchLength = math.min((b + 1)*valBatchSize - b*valBatchSize, valN - b*valBatchSize)

      local input = valInputs:narrow(1, batchStart, batchLength)
      local output = valOutputs:narrow(1, batchStart, batchLength)

      input = input:cuda()
      output = output:cuda()

      local valPreds = model:forward(input)
      local valCodes = mean.output
      local valLogVar = logvar.output

      local f = criterion:forward(valPreds, output)
      valProtocol[valIteration][2] = valProtocol[valIteration][2] + f
      valProtocol[valIteration][4] = valProtocol[valIteration][4] + KLD.loss

      valProtocol[valIteration][6] = valProtocol[valIteration][6] + torch.mean(torch.abs(valPreds - output))

      valProtocol[valIteration][8] = valProtocol[valIteration][8] + torch.mean(valCodes)
      valProtocol[valIteration][9] = valProtocol[valIteration][9] + torch.std(valCodes)
      valProtocol[valIteration][10] = valProtocol[valIteration][10] + torch.mean(valLogVar)

      local valPreds = decoder:forward(valCodes)

      local f = criterion:forward(valPreds, output)
      valProtocol[valIteration][3] = valProtocol[valIteration][3] + f
      valProtocol[valIteration][5] = valProtocol[valIteration][5] + KLD.loss

      valProtocol[valIteration][7] = valProtocol[valIteration][7] + torch.mean(torch.abs(valPreds - output))

      accValMeanPreds = lib.utils.appendTensor(accValMeanPreds, valPreds)
    end

    valProtocol[valIteration][1] = t
    for i = 2, valProtocol:size(2) do
      valProtocol[valIteration][i] = valProtocol[valIteration][i] / valNumBatches
    end

    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': validation loss ' .. valProtocol[valIteration][2])
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': KLD loss ' .. valProtocol[valIteration][4])
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': codes (mean) ' .. valProtocol[valIteration][8] .. ' / ' .. valProtocol[valIteration][9])
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': codes (logvar) ' .. valProtocol[valIteration][10])
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': validation error (abs) ' .. valProtocol[valIteration][6])

    local predFile = config['base_directory'] .. t .. '_predictions.h5'
    lib.utils.writeHDF5(predFile, accValMeanPreds)
    print('[Training] ' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': saved ' .. predFile)
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