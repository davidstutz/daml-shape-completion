require('torch')
require('nn')

--- @class PerChannelCriterion
local PerChannelCriterion, PerChannelCriterionParent = torch.class('nn.PerChannelCriterion', 'nn.Criterion')

--- Initialize.
function PerChannelCriterion:__init()
  self.criteria = nil
  self.weights = nil
end

--- Compute the loss.
-- @param input output of the network
-- @param target targets
-- @return loss
function PerChannelCriterion:updateOutput(input, target)
  assert(self.criteria)
  assert(self.weights)
  assert(#self.criteria == #self.weights)
  assert(input:dim() == target:dim())
  assert(input:size(2) == #self.criteria)
  assert(input:size(2) == target:size(2))

  self.output = 0
  for i = 1, #self.criteria do
    self.output = self.output + self.weights[i]*self.criteria[i]:forward(input:narrow(2, i, 1), target:narrow(2, i, 1))
  end

  return self.output
end

--- Compute the backward pass of the loss.
-- @param input original input to loss
-- @param target targets
-- @return gradients with respect to input
function PerChannelCriterion:updateGradInput(input, target)
  assert(self.criteria)
  assert(self.weights)
  assert(#self.criteria == #self.weights)
  assert(input:dim() == target:dim())
  assert(input:size(2) == #self.criteria)
  assert(input:size(2) == target:size(2))

  self.gradInput = input:clone():fill(0)
  for i = 1, #self.criteria do
    self.gradInput:narrow(2, i, 1):copy(self.weights[i]*self.criteria[i]:backward(input:narrow(2, i, 1), target:narrow(2, i, 1)))
  end

  return self.gradInput
end