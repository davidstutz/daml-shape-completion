require('torch')
require('nn')

--- @class GaussianToBernoulli
local GaussianToBernoulli, GaussianToBernoulliParent = torch.class('nn.GaussianToBernoulli', 'nn.Module')

--- Initialize.
function GaussianToBernoulli:__init(mean)
  self.logvar = -2
  self.eps = 0.5
end

--- Print dimensions of last layer.
-- @param input output of last layer
-- @return unchanged output of last layer
function GaussianToBernoulli:updateOutput(input)
  assert(self.logvar)

  -- for each pixel compute the probability of non-occupancy as
  -- \int_{-\infty}^\epsilon N(x|\mu,\sigma^2) dx
  -- i.e. as he cumulative density function of a Gaussian and
  -- a fixed parameter epsilon.
  -- So we compute
  -- p_i = 1 - \frac{1}{2}[1 + erf(\frac{\epsilon - \mu_i}{\sqrt{2} \sigma_i})]

  local function erf(x)
    -- constants
    local a1 =  0.254829592
    local a2 = -0.284496736
    local a3 =  1.421413741
    local a4 = -1.453152027
    local a5 =  1.061405429
    local p  =  0.3275911

    -- Save the sign of x
    local sign = torch.sign(x)
    x = torch.abs(x)

    -- A&S formula 7.1.26
    local ones = x:clone():fill(1)
    local t = torch.cdiv(ones, 1.0 + p*x)

    local y = torch.cmul(a5*t + a4, t)
    y = torch.cmul(y + a3, t)
    y = torch.cmul(y + a2, t)
    y = torch.cmul(y + a1, t)
    y = 1.0 - torch.cmul(y, torch.exp(torch.cmul(-x, x)))

    return torch.cmul(sign, y)
  end

  self.output = 0.5*(1 + erf(math.exp(-0.5*self.logvar)*0.707106781*(self.eps - input)))
  return self.output
end

--- Print the gradients of the next layer.
-- @param input original input of last layer
-- @param gradOutput gradients of next layer
-- @return unchanged gradients of next layer
function GaussianToBernoulli:updateGradInput(input, gradOutput)
  self.gradInput = torch.cmul(- math.exp(-0.5*self.logvar)/math.sqrt(2*math.pi)*torch.exp(-0.5*torch.exp(-self.logvar)*torch.cmul(self.eps - input, self.eps - input)), gradOutput)
  return self.gradInput
end