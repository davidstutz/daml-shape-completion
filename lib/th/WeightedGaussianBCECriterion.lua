require('torch')
require('nn')

--- @class WeightedGaussianBCECriterion
local WeightedGaussianBCECriterion, WeightedGaussianBCECriterionParent = torch.class('nn.WeightedGaussianBCECriterion', 'nn.Criterion')

--- Initialize.
function WeightedGaussianBCECriterion:__init()
  self.sizeAverage = false
  self.eps = 0
  self.weights = nil
  self.g2b = nil
  self.bce = nil
end

--- Compute the loss.
-- @param input output of the network
-- @param target targets
-- @return loss
function WeightedGaussianBCECriterion:updateOutput(input, target)
  assert(input:dim() == target:dim())
  assert(input:size(2) == target:size(2))
  assert(self.weights)

  if self.g2b == nil then
    self.g2b = nn.Sequential()
    local g2b = nn.GaussianToBernoulli()
    g2b.eps = self.eps
    self.g2b:add(g2b)
    self.g2b = self.g2b:cuda()
  end

  if self.bce == nil then
    self.bce = nn.WeightedBCECriterion()
    self.bce.sizeAverage = self.sizeAverage
    self.bce = self.bce:cuda()
  end

  self.bce.weights = self.weights
  self.bernoulli = self.g2b:forward(input)
  --require('image')
  --image.save('test.png', self.bernoulli[1][1])
  self.output = self.bce:forward(self.bernoulli, target)
  return self.output
end

--- Compute the backward pass of the loss.
-- @param input original input to loss
-- @param target targets
-- @return gradients with respect to input
function WeightedGaussianBCECriterion:updateGradInput(input, target)
  assert(input:dim() == target:dim())
  assert(input:size(2) == target:size(2))
  assert(self.weights)

  self.bce.weights = self.weights
  self.gradInput = self.bce:backward(self.bernoulli, target)
  self.gradInput = self.g2b:backward(input, self.gradInput)
  return self.gradInput
end