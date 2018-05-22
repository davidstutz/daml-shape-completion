-- Run model on config.json.

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

print('[Validation] reading ' .. configFile)
config = lib.utils.readJSON(configFile)

config['data_directory'] = config['data_directory'] .. '/'
if not lib.utils.directoryExists(config['data_directory']) then
  print('[Error] Data directory not found')
  os.exit()
end

if not lib.utils.directoryExists(config['base_directory']) then
  print('[Error] ' .. config['base_directory'] .. ' does not exist; train a model first - or, for KITTI, manually create the directory and place a trained model inside')
  os.exit()
end

-- Load data for validation.
valInputs_1 = checkAndRead(config['data_directory'] .. config['validation_inputs'])
valInputs_2 = checkAndRead(config['data_directory'] .. config['validation_lsdf_inputs'])
valInputs = torch.cat({valInputs_1:float(), valInputs_2}, 2)

-- ---------------------------------------------------------------
-- Read model.
-- ---------------------------------------------------------------

modelFile = config['base_directory'] .. config['inference_model_file']
if not lib.utils.fileExists(modelFile) then
  print('[Error] Model file ' .. modelFile .. ' not found.')
  print('[Error] Check the README to download the models and put the correct model in base_directory.')
  os.exit()
end

model = torch.load(modelFile)
decoder = model.modules[2]
print(model)

-- ---------------------------------------------------------------
-- Run model.
-- ---------------------------------------------------------------

local valN = valInputs:size(1)
local valBatchSize = 64
local valNumBatches = math.floor(valN/valBatchSize)

local accValPreds = nil
for b = 0, valNumBatches do
  local batchStart = b*valBatchSize + 1
  local batchLength = math.min((b + 1)*valBatchSize - b*valBatchSize, valN - b*valBatchSize)

  local input = valInputs:narrow(1, batchStart, batchLength):cuda()

  -- First pass to get mean, could be done faster but does the job.
  local valPreds = model:forward(input)
  local valCodes = model.modules[1].modules[#model.modules[1].modules - 2].modules[1].output

  -- Second pass with mean code.
  valPreds = decoder:forward(valCodes)
  accValPreds = lib.utils.appendTensor(accValPreds, valPreds)

  local time = os.date("*t")
  print('[Validation] ' .. (b + 1)*valBatchSize .. '|' .. valNumBatches*valBatchSize)
end

local predFile = config['base_directory'] .. '0_predictions.h5'
lib.utils.writeHDF5(predFile, accValPreds)
print('[Validation] ' .. predFile)