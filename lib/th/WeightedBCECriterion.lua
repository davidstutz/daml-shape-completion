require('torch')
require('nn')

--- @class WeightedBCECriterion
local WeightedBCECriterion, WeightedBCECriterionParent = torch.class('nn.WeightedBCECriterion', 'nn.Criterion')

function WeightedBCECriterion:__init()
  self.sizeAverage = false
  self.weights = nil
end

--- Compute the loss.
-- @param input output of the network
-- @param target targets
-- @return loss
function WeightedBCECriterion:updateOutput(input, target)
  assert(input:dim() == target:dim())
  assert(input:size(2) == target:size(2))
  assert(self.weights)

  self.output = - torch.cmul(target, torch.log(input + 1e-12)) - torch.cmul(1 - target, torch.log(1 - input + 1e-12))
  self.output = torch.sum(torch.cmul(self.weights, self.output))

  if self.sizeAverage then
    self.output = self.output / lib.utils.storageProd(#input)
  end

  return self.output
end

--- Compute the backward pass of the loss.
-- @param input original input to loss
-- @param target targets
-- @return gradients with respect to input
function WeightedBCECriterion:updateGradInput(input, target)
  assert(input:dim() == target:dim())
  assert(input:size(2) == target:size(2))
  assert(self.weights)

  self.gradInput = - torch.cdiv(target, input + 1e-12) + torch.cdiv(1 - target, 1 - input + 1e-12)
  self.gradInput = torch.cmul(self.weights, self.gradInput)

  if self.sizeAverage then
    self.gradInput = self.gradInput / lib.utils.storageProd(#input)
  end

  return self.gradInput
end
