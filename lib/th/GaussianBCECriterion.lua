require('torch')
require('nn')

--- @class GaussianBCECriterion
local GaussianBCECriterion, GaussianBCECriterionParent = torch.class('nn.GaussianBCECriterion', 'nn.Criterion')

--- Initialize.
function GaussianBCECriterion:__init()
  self.sizeAverage = false
  self.eps = 0
  self.g2b = nil
  self.bce = nil
end

--- Compute the loss.
-- @param input output of the network
-- @param target targets
-- @return loss
function GaussianBCECriterion:updateOutput(input, target)
  assert(input:dim() == target:dim())
  assert(input:size(2) == target:size(2))

  if self.g2b == nil then
    self.g2b = nn.Sequential()
    local g2b = nn.GaussianToBernoulli()
    g2b.eps = self.eps
    self.g2b:add(g2b)
    self.g2b = self.g2b:cuda()
  end

  if self.bce == nil then
    self.bce = nn.BCECriterion()
    self.bce.sizeAverage = self.sizeAverage
    self.bce = self.bce:cuda()
  end

  self.bernoulli = self.g2b:forward(input)
  self.output = self.bce:forward(self.bernoulli, target)
  return self.output
end

--- Compute the backward pass of the loss.
-- @param input original input to loss
-- @param target targets
-- @return gradients with respect to input
function GaussianBCECriterion:updateGradInput(input, target)
  assert(input:dim() == target:dim())
  assert(input:size(2) == target:size(2))

  self.gradInput = self.bce:backward(self.bernoulli, target)
  self.gradInput = self.g2b:backward(input, self.gradInput)
  return self.gradInput
end