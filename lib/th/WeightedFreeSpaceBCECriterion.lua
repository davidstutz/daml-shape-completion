require('torch')
require('nn')

--- @class WeightedFreeSpaceBCECriterion
local WeightedFreeSpaceBCECriterion, WeightedFreeSpaceBCECriterionParent = torch.class('nn.WeightedFreeSpaceBCECriterion', 'nn.Criterion')

--- Initialize.
function WeightedFreeSpaceBCECriterion:__init()
  self.sizeAverage = true
  self.criterion = nil
  self.weights = nil
end

--- Compute the loss.
-- @param input output of the network
-- @param target targets
-- @return loss
function WeightedFreeSpaceBCECriterion:updateOutput(input, target)
  assert(input:dim() == target:dim())
  assert(self.weights)

  if not self.criterion then
    self.criterion = nn.WeightedBCECriterion()
    self.criterion.weights = self.weights
    self.criterion.sizeAverage = self.sizeAverage
    self.criterion = self.criterion:cuda()
  end

  -- target equals 1 for free space
  self.criterion.weights = torch.cmul(self.weights:narrow(1, 1, target:size(1)), target)
  self.output = self.criterion:forward(input, 1 - target)

  return self.output
end

--- Compute the backward pass of the loss.
-- @param input original input to loss
-- @param target targets
-- @return gradients with respect to input
function WeightedFreeSpaceBCECriterion:updateGradInput(input, target)
  assert(input:dim() == target:dim())
  assert(self.criterion)

  self.criterion.weights = torch.cmul(self.weights, target)
  self.gradInput = self.criterion:backward(input, 1 - target)

  return self.gradInput
end