-- A more complex initialization module which allows non-uniform initialization not using the
-- layers' reset method.
require("nn")

--- Initialize a tensor with a fixed value.
-- @param tensor tensor to initialize
-- @param fanIn number of input units
-- @param fanOut number of output units
-- @param value value to initialize with
local function initFixed(tensor, fanIn, fanOut, value)
  if tensor then
    tensor:fill(value)
  end
end

--- Uniform initialization.
-- @param tensor tensor to initialize
-- @param fanIn number of input units
-- @param fanOut number of output units
-- @param value value to use as range for uniform intialization
local function initUniform(tensor, fanIn, fanOut, value)
  if tensor then
    tensor:uniform(-value, value)
  end
end

--- Initialize a tensor according to normal distribution.
-- @param tensor tensor to initialize
-- @param fanIn number of input units
-- @param fanOut number of output units
-- @param value value to initialize with
local function initNormal(tensor, fanIn, fanOut, value)
  if tensor then
    tensor:normal(0, value)
  end
end


--- Initialization scheme introduced in
-- "Efficient backprop", Yann Lecun, 1998
-- @param tensor tensor to initialize
-- @param fanIn number of input units
-- @param fanOut number of output units
-- @param value not used
local function initHeuristic(tensor, fanIn, fanOut, value)
  local std = math.sqrt(1/(3*fanIn))
  std = std * math.sqrt(3)

  if tensor then
    tensor:uniform(-std, std)
  end
end

--- Initialization scheme introduced in
-- "Understanding the difficulty of training deep feedforward neural networks", Xavier Glorot, 2010
-- @param tensor tensor to initialize
-- @param fanIn number of input units
-- @param fanOut number of output units
-- @param value not used
local function initXavier(tensor, fanIn, fanOut, value)
  local std = math.sqrt(2/(fanIn + fanOut))
  std = std * math.sqrt(3)

  if tensor then
    tensor:uniform(-std, std)
  end
end

--- Initialization scheme introduced in
-- "Understanding the difficulty of training deep feedforward neural networks", Xavier Glorot, 2010
-- @param tensor tensor to initialize
-- @param fanIn number of input units
-- @param fanOut number of output units
-- @param value not used
local function initXavierCaffe(tensor, fanIn, fanOut, value)
  local std = math.sqrt(1/fanIn)
  std = std * math.sqrt(3)

  if tensor then
    tensor:uniform(-std, std)
  end
end

--- Initialization scheme introduced in
-- @param tensor tensor to initialize
-- @param fanIn number of input units
-- @param fanOut number of output units
-- @param value not used
local function initKaiming(tensor, fanIn, fanOut, value)
  local std = math.sqrt(4/(fanIn + fanOut))
  std = std * math.sqrt(3)

  if tensor then
    tensor:uniform(-std, std)
  end
end

--- Get the init function by its name.
-- @param name name of the function
-- @return the function
local function getMethodByName(name)
  if name == 'fixed' then
    return initFixed
  elseif name == 'uniform' then
    return initUniform
  elseif name == 'normal' then
    return initNormal
  elseif name == 'heuristic' then
    return initHeuristic
  elseif name == 'xavier' then
    return initXavier
  elseif name == 'xavier_caffe' then
    return initXavierCaffe
  elseif name == 'kaiming' then
    return initKaiming
  else
    assert(false, 'unkown initialization method: ' .. name)
  end
end

local function initializeRecursive(model, weightsMethod, biasMethod, weightsValue, biasValue)
  assert(model)

  for i = 1, #model.modules do
    local m = model.modules[i]
    if m.modules ~= nil then
      initializeRecursive(m, weightsMethod, biasMethod, weightsValue, biasValue)
    elseif m.__typename == 'nn.SpatialConvolution' or
           m.__typename == 'nn.SpatialConvolutionMM' or
           m.__typename == 'nn.SpatialFullConvolution' then

      local fanIn = m.nInputPlane*m.kH*m.kW
      local fanOut = m.nOutputPlane*m.kH*m.kW
      weightsMethod(m.weight, fanIn, fanOut, weightsValue)
      biasMethod(m.bias, fanIn, fanOut, biasValue)
    elseif m.__typename == 'nn.VolumetricConvolution' or
           m.__typename == 'nn.VolumetricFullConvolution' or
           m.__typename == 'nn.VolumetricDilatedConvolution' then

      local fanIn = m.nInputPlane*m.kH*m.kW*m.kT
      local fanOut = m.nOutputPlane*m.kH*m.kW*m.kT
      weightsMethod(m.weight, fanIn, fanOut, weightsValue)
      biasMethod(m.bias, fanIn, fanOut, biasValue)
    elseif m.__typename == 'nn.Linear' then
      local fanIn = m.weight:size(2)
      local fanOut = m.weight:size(1)
      weightsMethod(m.weight, fanIn, fanOut, weightsValue)
      biasMethod(m.bias, fanIn, fanOut, biasValue)

    -- Explitly list layers that do no tneed to be initialized
    -- as to avoid forgetting to initialize newly createda layers.
    elseif m.__typename == 'nn.VolumetricBatchNormalization' or
           m.__typename == 'nn.VolumetricMaxPooling' or
           m.__typename == 'nn.VolumetricUpSamplingNearest' or
           m.__typename == 'nn.Sigmoid' or
           m.__typename == 'nn.ReLU' or
           m.__typename == 'nn.View' or
           m.__typename == 'nn.PrintDimensions' or
           m.__typename == 'nn.CheckNaN' or
           m.__typename == 'nn.GaussianKullbackLeiblerDivergence' or
           m.__typename == 'nn.GaussianReparameterizationSampler' or
           m.__typename == 'nn.Identity' or
           m.__typename == 'nn.PerChannelNonLinearity' or
           m.__typename == 'nn.MaximumLikelihoodPrior' then

      -- Nothing!
    else
      assert(false, 'weight initialization: unknown layer ' .. m.__typename) -- To be sure that all layers are initialized appropriately.
    end
  end
end

--- Use the given method to initialize all layers.
-- @param model model to initialize
-- @param methodName method to use
local function init(model, weightsMethodName, weightsValue, biasMethodName, biasValue)

  local weightsValue = weightsValue or 0.05
  local biasMethodName = biasMethodName or 'fixed'
  local biasValue = biasValue or 0.0

  local weightsMethod = getMethodByName(weightsMethodName)
  local biasMethod = getMethodByName(biasMethodName)

  -- loop over all convolutional modules
  initializeRecursive(model, weightsMethod, biasMethod, weightsValue, biasValue)
end

lib.init = init