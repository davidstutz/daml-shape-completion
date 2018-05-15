lib = {}

-- Utilities:
include('Utils.lua')

-- Layers and criteria:
include('PerChannelNonLinearity.lua')
include('PerChannelCriterion.lua')
include('GaussianKullbackLeiblerDivergence.lua')
include('GaussianReparameterizationSampler.lua')
include('WeightedBCECriterion.lua')
include('GaussianBCECriterion.lua')
include('WeightedGaussianBCECriterion.lua')
include('FreeSpaceBCECriterion.lua')
include('WeightedFreeSpaceBCECriterion.lua')
include('PointBCECriterion.lua')
include('GaussianFreeSpaceBCECriterion.lua')
include('WeightedGaussianFreeSpaceBCECriterion.lua')
include('GaussianPointBCECriterion.lua')
include('FixedVarianceGaussianNLLCriterion.lua')
include('MaximumLikelihoodPrior.lua')
include('MultipleCriterion.lua')
include('GaussianToBernoulli.lua')

-- Weight initialization:
include('WeightInitialization.lua')

-- Convolutional Models:
include('AutoEncoder.lua')
include('VariationalAutoEncoder.lua')

return lib