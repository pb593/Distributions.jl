struct AsymmetricLaplace <: ContinuousUnivariateDistribution
    μ::Real # location
    θ::Real # scale
    κ::Real # assymmetry

    AsymmetricLaplace(μ::Real, θ::Real, κ::Real) = new(μ, θ, κ)
end

mean(d::AsymmetricLaplace) = d.μ + (1 - d.κ^2) / (d.θ * d.κ)
var(d::AsymmetricLaplace) = (1 + d.κ^4) / (d.θ * d.κ)^2
skewness(d::AsymmetricLaplace) = 2 * (1 - d.κ^6) / (d.κ^4 + 1)^(3/2)
kurtosis(d::AsymmetricLaplace) = 6 * (1 + d.κ^8) / (1 +  d.κ^4)^2
median(d::AsymmetricLaplace) =
    d.μ +
    if (d.κ >= 1)
        (d.κ / d.θ) * log((1 + d.κ^2) / (2 * d.κ^2))
    else
        - log((1 + d.κ^2) / 2) / (d.θ * d.κ)
    end

entropy(d::AsymmetricLaplace) = 1 + log((1 + d.κ^2) / (d.κ * d.θ))
params(d::AsymmetricLaplace) = (d.μ, d.θ, d.κ)

rand(rng::AbstractRNG, d::AsymmetricLaplace) =
    rand(rng, Exponential(d.θ / d.κ)) - rand(rng, Exponential(d.θ * d.κ)) + d.μ
