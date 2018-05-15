require('torch')
require('nn')

--- @class PerChannelNonLinearity
local PerChannelNonLinearity, PerChannelNonLinearityParent = torch.class('nn.PerChannelNonLinearity', 'nn.Module')

--- Constructor.
function PerChannelNonLinearity:__init()
  self.layers = nil
end

--- Forward pass.
-- @param input input to layer
-- @return output of layer
function PerChannelNonLinearity:updateOutput(input)
  assert(self.layers)
  assert(#self.layers == input:size(2))

  self.output = input:clone()
  for i = 1, #self.layers do
    self.output:narrow(2, i, 1):copy(self.layers[i]:forward(input:narrow(2, i, 1)))
  end

  return self.output
end

--- Backward pass.
-- @param input input to layer
-- @param gradOutput gradients from top layer
-- @return gradients from this layer
function PerChannelNonLinearity:updateGradInput(input, gradOutput)
  assert(self.layers)
  assert(#self.layers == input:size(2))
  assert(gradOutput:size(2) == input:size(2))

  self.gradInput = input:clone():fill(0)
  for i = 1, #self.layers do
    self.gradInput:narrow(2, i, 1):copy(self.layers[i]:backward(input:narrow(2, i, 1), gradOutput:narrow(2, i, 1)))
  end

  return self.gradInput
end