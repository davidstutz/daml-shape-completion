require('torch')
require('nn')

--- @class PointBCECriterion
local PointBCECriterion, PointBCECriterionParent = torch.class('nn.PointBCECriterion', 'nn.Criterion')

--- Initialize.
function PointBCECriterion:__init()
  self.sizeAverage = true
  self.criterion = nil
end

--- Compute the loss.
-- @param input output of the network
-- @param target targets
-- @return loss
function PointBCECriterion:updateOutput(input, target)
  assert(input:dim() == target:dim())

  if not self.criterion then
    self.criterion = nn.WeightedBCECriterion()
    self.criterion.sizeAverage = self.sizeAverage
    self.criterion = self.criterion:cuda()
  end

  self.criterion.weights = target
  self.output = self.criterion:forward(input, target)

  return self.output
end

--- Compute the backward pass of the loss.
-- @param input original input to loss
-- @param target targets
-- @return gradients with respect to input
function PointBCECriterion:updateGradInput(input, target)
  assert(input:dim() == target:dim())
  assert(self.criterion)

  self.criterion.weights = target
  self.gradInput = self.criterion:backward(input, target)

  return self.gradInput
end