require('torch')
require('nn')

--- @class GaussianReparameterizationSampler
local GaussianReparameterizationSampler, GaussianReparameterizationSamplerParent = torch.class('nn.GaussianReparameterizationSampler', 'nn.Module')

--- Initialize.
function GaussianReparameterizationSampler:__init()
  -- Nothing ...
end

--- Sample from the provided mean and variance using the reparameterization trick.
-- @param input table of two elements, mean and log variance
-- @return sample
function GaussianReparameterizationSampler:updateOutput(input)
  assert(#input == 2)

  local mean, logVar = table.unpack(input)
  self.eps = torch.randn(input[1]:size()):cuda()
  self.output = torch.cmul(torch.exp(0.5*logVar), self.eps) + mean

  return self.output
end

--- Backward pass of the sampler.
-- @param input table of two elements, mean and log variance
-- @param gradOutput gradients of top layer
-- @return gradients with respect to input, table of two elements
function GaussianReparameterizationSampler:updateGradInput(input, gradOutput)
  self.gradInput = {}

  local mean, logVar = table.unpack(input)
  self.gradInput[1] = gradOutput
  self.gradInput[2] = torch.cmul(torch.cmul(0.5*torch.exp(0.5*logVar), self.eps), gradOutput)

  return self.gradInput
end