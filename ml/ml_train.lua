-- Maximum liklihood for inference

require('torch')
require('cutorch')
require('nn')
require('nnx')
require('cunn')
require('cunnx')
require('optim')
require('hdf5')
require('lfs')
require('image')

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

-- For later simplicity.
N = valInputs:size(1)
channels = valInputs:size(2)
height = valInputs:size(3)
width = valInputs:size(4)
depth = valInputs:size(5)

model = torch.load(config['base_directory'] .. config['prior_model_file'])

-- Get the decoder of the pre-trained model.
MLP = nn.MaximumLikelihoodPrior()
MLP.sizeAverage = false
MLP.lambda = config['prior_weight']

decoder = nn.Sequential()
decoder:add(MLP)
decoder:add(model.modules[2]) -- !
print(decoder)

codeSize = model.modules[2].modules[1].weight:size(2)
batchSize = config['batch_size']



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

-- Learning hyperparameters.
iterations = config['iterations']
minimumLearningRate = config['minimum_learning_rate']
learningRateDecay = config['decay_learning_rate']
maximumMomentum = config['maximum_momentum']
momentumDecay = config['decay_momentum']
decayIterations = config['decay_iterations']
iterations = config['iterations']
lossIterations = config['loss_iterations']

-- https://github.com/torch/cutorch/issues/497
predictions = torch.Tensor(N, batchSize, channels, height, width, depth)
protocol = torch.Tensor(iterations, N, 5):fill(0)

for n = 1, 5 do
  learningRate = config['learning_rate']
  momentum = config['momentum']

  -- Update state with learning rate and momentum.
  state = {
    learningRate = learningRate,
    momentum = momentum,
    learningRateDecay = 0 -- will be done manually below
  }

  print('[Training] ' .. n .. ': reset learning rate ' .. state.learningRate)
  print('[Training] ' .. n .. ': reset momentum ' .. state.momentum)

  local code = torch.randn(batchSize, codeSize)
  code:fill(0)
  code = code:cuda()

  local gradCode = torch.Tensor(batchSize, codeSize):fill(0)
  gradCode = gradCode:cuda()

  local points = torch.Tensor(batchSize, 1, height, width, depth)
  local space = torch.Tensor(batchSize, 1, height, width, depth)
  local output = nil

  if valOutputs then
    output = torch.Tensor(batchSize, 2, height, width, depth)
  end

  for i = 1, batchSize do
    points[i]:copy(valPoints[n])
    space[i]:copy(valSpace[n])

    if output then
      output[i]:copy(valOutputs[n])
    end
  end

  points = points:cuda()
  space = space:cuda()

  if output then
    output = output:cuda()
  end

  for t = 1, iterations do
    protocol[t][n][1] = t

    --- Definition of the objective on the current mini-batch.
    -- This will be the objective fed to the optimization algorithm.
    -- @param x input parameters
    -- @return object value, gradients
    local feval = function(x)

      if x ~= code then
        code:copy(x)
      end

      local pred = decoder:forward(code)
      local f = criterion:forward(pred, {points, space, points, space})

      local df_do = criterion:backward(pred, {points, space, points, space})
      gradCode = decoder:backward(code, df_do)

      protocol[t][n][2] = f
      protocol[t][n][3] = MLP.loss

      if output then
        protocol[t][n][5] = torch.mean(torch.abs(pred - output))
      end

      f = f + MLP.loss

      return f, gradCode
    end

    _, _ = optim.sgd(feval, code, state)
    local time = os.date("*t")

    if t%lossIterations == 0 then
      local smoothedLoss = torch.mean(protocol:narrow(1, t - lossIterations + 1, lossIterations):narrow(2, n, 1):narrow(3, 2, 1))
      local smoothedKLD = torch.mean(protocol:narrow(1, t - lossIterations + 1, lossIterations):narrow(2, n, 1):narrow(3, 3, 1))
      local smoothedError = torch.mean(protocol:narrow(1, t - lossIterations + 1, lossIterations):narrow(2, n, 1):narrow(3, 5, 1))
      local loss = smoothedLoss + smoothedKLD

      if t > 2*lossIterations then
        local smoothedLastLoss = torch.mean(protocol:narrow(1, t - 2*lossIterations + 1, lossIterations):narrow(2, n, 1):narrow(3, 2, 1))
        local smoothedLastKLD = torch.mean(protocol:narrow(1, t - 2*lossIterations + 1, lossIterations):narrow(2, n, 1):narrow(3, 3, 1))
        local lastLoss = smoothedLastLoss + smoothedLastKLD

        if math.abs(loss -lastLoss) < 0.001 then
          print('[Training] ' .. n .. '|' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': breaking')
          break
        end
      end

      print('[Training] ' .. n .. '|' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec ..': ' .. loss .. ' / ' .. smoothedLoss .. ' / ' .. smoothedKLD .. ' / ' .. smoothedError)
    end

    -- Decay learning rate and KLD weight, do this before resetting all smoothed
    -- statistics.
    if t%decayIterations == 0 then
      state.learningRate = math.max(minimumLearningRate, state.learningRate*learningRateDecay)
      state.momentum = math.min(maximumMomentum, state.momentum*momentumDecay)
      print('[Training] ' .. n .. '|' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': learning rate ' .. state.learningRate)
      print('[Training] ' .. n .. '|' .. t .. '|' .. iterations .. '|' .. time.hour .. ':' .. time.min .. ':' .. time.sec .. ': momentum ' .. state.momentum)
    end

    -- !
    collectgarbage()
  end

  predictions[n] = decoder:forward(code):float()
end

lib.utils.writeHDF5(config['base_directory'] .. '0_predictions.h5', predictions)
print('[Training] wrote ' .. config['base_directory'] .. '0_predictions.h5')
lib.utils.writeHDF5(config['base_directory'] .. config['protocol_file'], protocol)
print('[Training] protocol ' .. config['base_directory'] .. config['protocol_file'])