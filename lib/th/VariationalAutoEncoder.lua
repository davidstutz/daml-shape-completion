-- Implementation of simple convolutional encoder/decoder achitecture with
-- variable number of channels, layers and kernel sizes.

require('nn')
require('cunn')
require('nnx')
require('cunnx')

local models = {}

--- Default options for the auto-encoder, encoder and decoder models.
models.config = {
  height = 0, -- input height
  width = 0, -- input width
  depth = 0, -- input depth
  code = 0, -- code dimension to use
  channels = {}, -- number of convolutional layers and number of channels for each
  kernelSizes = {}, -- kernel sizes of the convolutional layers indicated by channels
  pooling = nil, -- where to put pooling
  transfer = nn.ReLU, -- a function that only takes one argument, which is "inplace"
  transfers = nil, -- where to put transfer functions, indicated by true/false
  normalizations = nil, -- where to put batch normalization
  dataChannels = 1, -- number of input channels
  outputChannels = 1, -- number of output channels
  poolingSizes = nil, -- pooling kernel sizes
  printDimensions = false, -- whether to print dimensions after each layer
  checkNaN = false, -- whether to check for NaN values after each layer
}

--- Simple encoder structure as also explained by models.autoEncoder.
-- @param model model to add encoder to
-- @param config configuration as illustrated in models.autoEncoderConfig
-- @return model
function models.encoder(model, config)
  local context = {}
  model, context = lib.autoEncoder.encoder(model, config)
  model:remove(#model.modules) -- Remove last linear layer to replace by mean and logvar

  context.mean = nn.Linear(context.hidden, config.code)
  context.logVar = nn.Linear(context.hidden, config.code)

  local meanLogVar = nn.ConcatTable()
  meanLogVar:add(context.mean)
  meanLogVar:add(context.logVar)
  model:add(meanLogVar)

  -- KLD is useful for setting weights and getting loss!
  context.KLD = nn.GaussianKullbackLeiblerDivergence()
  model:add(context.KLD)
  context.sampler = nn.GaussianReparameterizationSampler()
  model:add(context.sampler)

  return model, context
end

--- Simple decoder structure as also explained by models.autoEncoder.
-- @param model model to add decoder to
-- @param config configuration as illustrated in models.autoEncoderConfig
-- @return model
function models.decoder(model, config)
  local model, context = lib.autoEncoder.decoder(model, config)
  return model, context
end

--- Sets up a decoder/encoder architecture with the given code dimensionality,
-- number of channels for each layer and the corresponding kernel sizes.
-- @param model model to add encoder and decoder to
-- @param config configuration as illustrated in models.autoEncoderConfig
-- @return model
function models.autoEncoder(model, config)
  local model = model or nn.Sequential()

  local context = {}
  local encoder = nn.Sequential()
  encoder, context = models.encoder(encoder, config)

  local decoder = nn.Sequential()
  decoder, _ = models.decoder(decoder, config)

  model:add(encoder)
  model:add(decoder)

  context.encoder = encoder
  context.decoder = decoder
  return model, context
end

lib.variationalAutoEncoder = models