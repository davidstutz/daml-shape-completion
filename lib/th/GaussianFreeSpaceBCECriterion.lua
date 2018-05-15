require('torch')
require('nn')

--- @class GaussianFreeSpaceBCECriterion
local GaussianFreeSpaceBCECriterion, GaussianFreeSpaceBCECriterionParent = torch.class('nn.GaussianFreeSpaceBCECriterion', 'nn.Criterion')

--- Initialize.
function GaussianFreeSpaceBCECriterion:__init()
  self.sizeAverage = true
  self.criterion = nil
end

--- Compute the loss.
-- @param input output of the network
-- @param target targets
-- @return loss
function GaussianFreeSpaceBCECriterion:updateOutput(input, target)
  assert(input:dim() == target:dim())

  if not self.criterion then
    self.criterion = nn.WeightedGaussianBCECriterion()
    self.criterion.sizeAverage = self.sizeAverage
    self.criterion = self.criterion:cuda()
  end

  self.criterion.weights = target
  self.output = self.criterion:forward(input, 1 - target)

  return self.output
end

--- Compute the backward pass of the loss.
-- @param input original input to loss
-- @param target targets
-- @return gradients with respect to input
function GaussianFreeSpaceBCECriterion:updateGradInput(input, target)
  assert(input:dim() == target:dim())
  assert(self.criterion)

  self.criterion.weights = target
  self.gradInput = self.criterion:backward(input, 1 - target)

  return self.gradInput
end