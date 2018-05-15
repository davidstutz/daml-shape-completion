require('torch')
require('nn')

--- @class FixedVarianceGaussianNLLCriterion
local FixedVarianceGaussianNLLCriterion, FixedVarianceGaussianNLLCriterionParent = torch.class('nn.FixedVarianceGaussianNLLCriterion', 'nn.Criterion')

--- Initialize.
function FixedVarianceGaussianNLLCriterion:__init()
  self.logvar = -2
  self.sizeAverage = true
end

--- Compute the loss.
-- @param input output of the network
-- @param target targets
-- @return loss
function FixedVarianceGaussianNLLCriterion:updateOutput(input, target)
  assert(input:dim() == target:dim())
  assert(input:size(2) == target:size(2))

  self.output = 0.5*math.exp(-self.logvar)*torch.cmul(target - input, target - input)
  self.output = torch.sum(self.output)

  if self.sizeAverage then
    self.output = self.output / lib.utils.storageProd(#input)
  end

  return self.output
end

--- Compute the backward pass of the loss.
-- @param input original input to loss
-- @param target targets
-- @return gradients with respect to input
function FixedVarianceGaussianNLLCriterion:updateGradInput(input, target)

  self.gradInput = -1*(target - input)*torch.exp(-self.logvar)

  if self.sizeAverage then
    self.gradInput = self.gradInput / lib.utils.storageProd(#input)
  end

  return self.gradInput
end