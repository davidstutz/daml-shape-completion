require('nn')
require('torch')

--- @class MaximumLikelihoodPrior
local MaximumLikelihoodPrior, MaximumLikelihoodPriorParent = torch.class('nn.MaximumLikelihoodPrior', 'nn.Module')

--- Initialize.
-- @param lambda weight of loss
function MaximumLikelihoodPrior:__init(lambda)
  self.lambda = lambda or 1
  self.sizeAverage = false
  self.loss = nil
end

--- Compute the Kullback-Leiber divergence; however, the input remains
-- unchanged - the divergence is saved in KullBackLeiblerDivergence.loss.
-- @param input table of two elements, mean and log variance
-- @param table of wo elements, mean and log variance
function MaximumLikelihoodPrior:updateOutput(input)
  -- + math.log(math.pow(2*math.pi, input:size(1)/2.))
  self.loss = self.lambda * (math.log(math.pow(2*math.pi, input:size(1)/2.)) + 0.5*torch.sum(torch.cmul(input, input)))

  if self.sizeAverage then
    self.loss = self.loss/lib.utils.storageProd(#input)
  end

  self.output = input
  return self.output
end

--- Compute the backward pass of the Kullback-Leibler Divergence.
-- @param input original inpur as table of two elements, mean and log variance
-- @param gradOutput gradients from top layer, table of two elements, mean and log variance
-- @param gradients with respect to input, table of two elements
function MaximumLikelihoodPrior:updateGradInput(input, gradOutput)
  local gradInput = self.lambda * input

  if self.sizeAverage then
    gradInput = gradInput/lib.utils.storageProd(#input)
  end

  self.gradInput = gradInput + gradOutput
  return self.gradInput
end