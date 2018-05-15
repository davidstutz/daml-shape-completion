require('torch')
require('nn')

--- @class WeightedGaussianFreeSpaceBCECriterion
local WeightedGaussianFreeSpaceBCECriterion, WeightedGaussianFreeSpaceBCECriterionParent = torch.class('nn.WeightedGaussianFreeSpaceBCECriterion', 'nn.Criterion')

--- Initialize.
function WeightedGaussianFreeSpaceBCECriterion:__init()
  self.sizeAverage = true
  self.criterion = nil
  self.weights = nil
end

--- Compute the loss.
-- @param input output of the network
-- @param target targets
-- @return loss
function WeightedGaussianFreeSpaceBCECriterion:updateOutput(input, target)
  assert(input:dim() == target:dim())

  if not self.criterion then
    self.criterion = nn.WeightedGaussianBCECriterion()
    self.criterion.sizeAverage = self.sizeAverage
    self.criterion = self.criterion:cuda()
  end

  self.criterion.weights = torch.cmul(self.weights:narrow(1, 1, target:size(1)), target)
  self.output = self.criterion:forward(input, 1 - target)

  return self.output
end

--- Compute the backward pass of the loss.
-- @param input original input to loss
-- @param target targets
-- @return gradients with respect to input
function WeightedGaussianFreeSpaceBCECriterion:updateGradInput(input, target)
  assert(input:dim() == target:dim())
  assert(self.criterion)

  self.criterion.weights = torch.cmul(self.weights, target)
  self.gradInput = self.criterion:backward(input, 1 - target)

  return self.gradInput
end