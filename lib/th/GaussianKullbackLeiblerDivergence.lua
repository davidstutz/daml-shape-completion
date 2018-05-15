require('nn')
require('torch')

--- @class GaussianKullbackLeiblerDivergence
local GaussianKullbackLeiblerDivergence, GaussianKullbackLeiblerDivergenceParent = torch.class('nn.GaussianKullbackLeiblerDivergence', 'nn.Module')

--- Initialize.
-- @param lambda weight of loss
function GaussianKullbackLeiblerDivergence:__init(lambda)
  self.lambda = lambda or 1
  self.lambda2 = 1
  self.sizeAverage = true
  self.loss = 0
end

--- Compute the Kullback-Leiber divergence; however, the input remains
-- unchanged - the divergence is saved in KullBackLeiblerDivergence.loss.
-- @param input table of two elements, mean and log variance
-- @param table of wo elements, mean and log variance
function GaussianKullbackLeiblerDivergence:updateOutput(input)
  assert(#input == 2)

  -- Save the loss for monitoring.
  local mean, logVar = table.unpack(input)
  self.loss = self.lambda * 0.5 * torch.sum(torch.pow(mean, 2) + self.lambda2*(torch.exp(logVar) - 1 - logVar))

  if self.sizeAverage then
    self.loss = self.loss/(lib.utils.storageProd(#input[1]))
  end

  self.output = input
  return self.output
end

--- Compute the backward pass of the Kullback-Leibler Divergence.
-- @param input original inpur as table of two elements, mean and log variance
-- @param gradOutput gradients from top layer, table of two elements, mean and log variance
-- @param gradients with respect to input, table of two elements
function GaussianKullbackLeiblerDivergence:updateGradInput(input, gradOutput)
  assert(#gradOutput == 2)

  local mean, logVar = table.unpack(input)
  self.gradInput = {}
  self.gradInput[1] = self.lambda*mean
  self.gradInput[2] = self.lambda*self.lambda2*0.5*(torch.exp(logVar) - 1)

  if self.sizeAverage then
    self.gradInput[1] = self.gradInput[1]/(lib.utils.storageProd(#input[1]))
    self.gradInput[2] = self.gradInput[2]/(lib.utils.storageProd(#input[2]))
  end

  self.gradInput[1] = self.gradInput[1] + gradOutput[1]
  self.gradInput[2] = self.gradInput[2] + gradOutput[2]

  return self.gradInput
end