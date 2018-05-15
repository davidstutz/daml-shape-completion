require('torch')
require('nn')

--- @class MultipleCriterion
local MultipleCriterion, MultipleCriterionParent = torch.class('nn.MultipleCriterion', 'nn.Criterion')

--- Initialize.
function MultipleCriterion:__init()
  self.criteria = nil
  self.weights = nil
  self.channels = nil
end

--- Compute the loss.
-- @param input output of the network
-- @param target targets
-- @return loss
function MultipleCriterion:updateOutput(input, target)
  assert(self.criteria)
  assert(self.weights)
  assert(#self.criteria > 0)
  assert(#self.weights == #self.criteria)
  assert(#target == #self.criteria)

  self.output = 0
  for i = 1, #self.criteria do
    local input_i = input
    if self.channels ~= nil then
      input_i = input:narrow(2, self.channels[i], 1)
    end

    self.output = self.output + self.weights[i]*self.criteria[i]:forward(input_i, target[i])
  end

  return self.output
end

--- Compute the backward pass of the loss.
-- @param input original input to loss
-- @param target targets
-- @return gradients with respect to input
function MultipleCriterion:updateGradInput(input, target)
  assert(self.criteria)
  assert(self.weights)
  assert(#self.criteria > 0)
  assert(#self.weights == #self.criteria)
  assert(#target == #self.criteria)

  self.gradInput = input:clone():fill(0)
  for i = 1, #self.criteria do
    local input_i = input
    if self.channels ~= nil then
      input_i = input:narrow(2, self.channels[i], 1)
    end

    if self.channels ~= nil then
      self.gradInput:narrow(2, self.channels[i], 1):copy(self.gradInput:narrow(2, self.channels[i], 1) + self.weights[i]*self.criteria[i]:backward(input_i, target[i]))
    else
      self.gradInput = self.gradInput + self.weights[i]*self.criteria[i]:backward(input_i, target[i])
    end
  end

  return self.gradInput
end