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

function models.computeHeightsWidthsDepths(config)
  assert(config.height > 0)
  assert(config.width > 0)
  assert(config.depth > 0)
  assert(config.poolingSizes == nil or #config.poolingSizes == #config.channels)
  assert(config.pooling == nil or #config.pooling == #config.channels)
  assert(config.strides == nil or #config.strides == #config.channels)

  local channels = config.channels
  local pooling = config.pooling
  local kernelSizes = config.kernelSizes
  local padding = config.padding
  local height = config.height
  local width = config.width
  local depth = config.depth
  local poolingSizes = config.poolingSizes
  local heights = {}
  local widths = {}
  local depths = {}

  for i = 1, #channels do
    heights[i] = height
    widths[i] = width
    depths[i] = depth

    if i > 1 then
      heights[i] = heights[i - 1]
      widths[i] = widths[i - 1]
      depths[i] = depths[i - 1]
    end

    if padding then
      local paddingSize = math.floor(kernelSizes[i]/2)
      heights[i] = heights[i] - 2*paddingSize
      widths[i] = widths[i] - 2*paddingSize
      depths[i] = depths[i] - 2*paddingSize
    end

    if pooling == nil or pooling[i] then
      local poolingSize = {2, 2, 2 }

      if poolingSizes ~= nil then
        if type(poolingSizes[i]) ~= 'table' then
          poolingSize = {poolingSizes[i], poolingSizes[i], poolingSizes[i] }
        else
          poolingSize = poolingSizes[i]
        end
      end

      heights[i] = math.floor(heights[i]/poolingSize[1])
      widths[i] = math.floor(widths[i]/poolingSize[2])
      depths[i] = math.floor(depths[i]/poolingSize[3])
    end

    if strides ~= nil and strides[i] > 1 then
      assert(strides[i] <= 2)
      heights[i] = math.floor(heights[i]/strides[i])
      widths[i] = math.floor(widths[i]/strides[i])
      depths[i] = math.floor(depths[i]/strides[i])
    end
  end

  return heights, widths, depths
end

--- Simple encoder structure as also explained by models.autoEncoder.
-- @param model model to add encoder to
-- @param config configuration as illustrated in models.autoEncoderConfig
-- @return model
function models.encoder(model, config)
  assert(config.height > 0)
  assert(config.width > 0)
  assert(config.depth > 0)
  assert(config.code > 0)
  assert(#config.channels > 0)
  assert(#config.channels == #config.kernelSizes)
  assert(config.poolingSizes == nil or #config.poolingSizes == #config.channels)
  assert(config.pooling == nil or #config.pooling == #config.channels)
  assert(config.transfers == nil or #config.channels == #config.transfers)
  assert(config.transfer)
  assert(config.normalizations == nil or #config.normalizations == #config.channels)
  assert(config.dataChannels > 0)

  local code = config.code
  local channels = config.channels
  local kernelSizes = config.kernelSizes
  local pooling = config.pooling
  local poolingSizes = config.poolingSizes
  local transfers = config.transfers
  local transfer = config.transfer
  local normalizations = config.normalizations
  local dataChannels = config.dataChannels
  local printDimensions = config.printDimensions
  local checkNaN = config.checkNaN

  -- Keeping track of the image/feature map sizes.
  local context = {}
  context.heights, context.widths, context.depths = models.computeHeightsWidthsDepths(config)

  -- Encoder part.
  for i = 1, #channels do

    -- Determine input channels.
    local inChannels = dataChannels
    if i > 1 then
      inChannels = channels[i - 1]
    end

    -- Output channels and kernels.
    local outChannels = channels[i]
    local kernelSize = kernelSizes[i]
    local padding = math.floor(kernelSize/2)

    local convolution = nn.VolumetricConvolution(inChannels, outChannels, kernelSize, kernelSize, kernelSize, 1, 1, 1, padding, padding, padding)

    if normalizations ~= nil and normalizations[i] then
      convolution:noBias()
    end

    model:add(convolution)

    if printDimensions then model:add(nn.PrintDimensions()) end
    if checkNaN then model:add(nn.CheckNaN()) end

    if normalizations ~= nil and normalizations[i] then
      model:add(nn.VolumetricBatchNormalization(outChannels))
      if printDimensions then model:add(nn.PrintDimensions()) end
      if checkNaN then model:add(nn.CheckNaN()) end
    end

    if transfers == nil or transfers[i] then
      model:add(transfer(true))
      if printDimensions then model:add(nn.PrintDimensions()) end
      if checkNaN then model:add(nn.CheckNaN()) end
    end

    if pooling == nil or pooling[i] then
      local poolingSize =  {2, 2, 2}

      if poolingSizes ~= nil then
        if type(poolingSizes[i]) ~= 'table' then
          poolingSizes[i] = {poolingSizes[i], poolingSizes[i], poolingSizes[i]}
        else
          poolingSize = poolingSizes[i]
        end
      end

      model:add(nn.VolumetricMaxPooling(poolingSize[1], poolingSize[3], poolingSize[2], poolingSize[1], poolingSize[3], poolingSize[2], 0, 0, 0))

      if printDimensions then model:add(nn.PrintDimensions()) end
      if checkNaN then model:add(nn.CheckNaN()) end
    end
  end

  local codeChannels = channels[#channels]
  local codeHeight = context.heights[#context.heights]
  local codeWidth = context.widths[#context.widths]
  local codeDepth = context.depths[#context.depths]

  context.hidden = codeChannels*codeHeight*codeWidth*codeDepth
  assert(context.hidden > 0)

  model:add(nn.View(context.hidden))
  if printDimensions then model:add(nn.PrintDimensions()) end
  if checkNaN then model:add(nn.CheckNaN()) end

  context.code = nn.Linear(context.hidden, code)
  model:add(context.code)
  -- Not after last layer!

  return model, context
end

--- Simple decoder structure as also explained by models.autoEncoder.
-- @param model model to add decoder to
-- @param config configuration as illustrated in models.autoEncoderConfig
-- @return model
function models.decoder(model, config)
  assert(config.height > 0)
  assert(config.width > 0)
  assert(config.depth > 0)
  assert(config.code > 0)
  assert(#config.channels > 0)
  assert(config.outputChannels > 0)
  assert(#config.channels == #config.kernelSizes)
  assert(config.poolingSizes == nil or #config.poolingSizes == #config.channels)
  assert(config.pooling == nil or #config.pooling == #config.channels)
  assert(config.transfers == nil or #config.channels == #config.transfers)
  assert(config.transfer)
  assert(config.normalizations == nil or #config.normalizations == #config.channels)
  assert(config.strides == nil or #config.strides == #config.channels)
  assert(config.dilations == nil or #config.dilations == #config.channels)

  local code = config.code
  local channels = config.channels
  local outputChannels = config.outputChannels
  local kernelSizes = config.kernelSizes
  local pooling = config.pooling
  local poolingSizes = config.poolingSizes
  local padding = config.padding
  local transfers = config.transfers
  local transfer = config.transfer
  local normalizations = config.normalizations
  local printDimensions = config.printDimensions
  local checkNaN = config.checkNaN

  -- Determine heights and widths.
  local context = {}
  context.heights, context.widths, context.depths = models.computeHeightsWidthsDepths(config)

  local codeChannels = channels[#channels]
  local codeHeight = context.heights[#context.heights]
  local codeWidth = context.widths[#context.widths]
  local codeDepth = context.depths[#context.depths]

  context.hidden = codeChannels*codeHeight*codeWidth*codeDepth
  model:add(nn.Linear(config.code, context.hidden))
  if printDimensions then model:add(nn.PrintDimensions()) end
  if checkNaN then model:add(nn.CheckNaN()) end

  model:add(nn.View(codeChannels, codeHeight, codeWidth, codeDepth))
  if printDimensions then model:add(nn.PrintDimensions()) end
  if checkNaN then model:add(nn.CheckNaN()) end

  -- Decoder part.
  for i = #channels, 1, -1 do

    if pooling == nil or pooling[i] then
      local poolingSize = {2, 2, 2}

      if poolingSizes ~= nil then
        if type(poolingSizes[i]) ~= 'table' then
          poolingSizes[i] = {poolingSizes[i], poolingSizes[i], poolingSizes[i] }
        else
          poolingSize = poolingSizes[i]
        end
      end

      model:add(nn.VolumetricUpSamplingNearest(poolingSize[1], poolingSize[2], poolingSize[3]))
      if printDimensions then model:add(nn.PrintDimensions()) end
      if checkNaN then model:add(nn.CheckNaN()) end
    end

    -- Determine input channels.
    local outChannels = outputChannels
    if i > 1 then
      outChannels = channels[i - 1]
    end

    -- Output channels and kernels.
    local inChannels = channels[i]
    local kernelSize = kernelSizes[i]
    local padding = math.floor(kernelSize/2)

    local convolution = nn.VolumetricConvolution(inChannels, outChannels, kernelSize, kernelSize, kernelSize, 1, 1, 1, padding, padding, padding)

    if normalizations ~= nil and normalizations[i] and i > 1 then
      convolution:noBias()
    end

    model:add(convolution)

    if printDimensions then model:add(nn.PrintDimensions()) end
    if checkNaN then model:add(nn.CheckNaN()) end

    -- not sure whether i > 1 should also be checked
    if normalizations ~= nil and normalizations[i] and i > 1 then
      model:add(nn.VolumetricBatchNormalization(outChannels))
      if printDimensions then model:add(nn.PrintDimensions()) end
      if checkNaN then model:add(nn.CheckNaN()) end
    end

    if transfers == nil or transfers[i] or i == 1 then
      if i > 1 then
        model:add(transfer(true))
        if printDimensions then model:add(nn.PrintDimensions()) end
        if checkNaN then model:add(nn.CheckNaN()) end
      else
        if centerMean ~= nil then
          context.uncenter = nn.UnCenter(centerMean)
          model:add(context.uncenter)
        end

        --model:add(nn.Sigmoid(true))
        -- Not after last layer!
      end
    end
  end

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

  context['encoder'] = encoder
  context['decoder'] = decoder
  return model, context
end

lib.autoEncoder = models