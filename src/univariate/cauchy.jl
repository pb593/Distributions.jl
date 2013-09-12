immutable Cauchy <: ContinuousUnivariateDistribution
    location::Float64
    scale::Float64
    function Cauchy(l::Real, s::Real)
        s > zero(s) || error("scale must be positive")
	new(float64(l), float64(s))
    end
end

Cauchy(l::Real) = Cauchy(l, 1.0)
Cauchy() = Cauchy(0.0, 1.0)

## Support
@continuous_distr_support Cauchy -Inf Inf

## Properties
mean(d::Cauchy) = NaN
median(d::Cauchy) = d.location
mode(d::Cauchy) = d.location
modes(d::Cauchy) = [mode(d)]

var(d::Cauchy) = NaN
skewness(d::Cauchy) = NaN
kurtosis(d::Cauchy) = NaN

entropy(d::Cauchy) = log(d.scale) + log(4.0 * pi)

## Functions
pdf(d::Cauchy, x::Real) = 1/(pi*d.scale*(1+((x-d.location)/d.scale)^2))
logpdf(d::Cauchy, x::Real) = -log(pi) - log(d.scale) - log1p(((x-d.location)/d.scale)^2)

cdf(d::Cauchy, x::Real) = atan2(one(x),-(x-d.location)/d.scale)/pi
ccdf(d::Cauchy, x::Real) = atan2(one(x),(x-d.location)/d.scale)/pi

quantile(d::Cauchy, p::Real) = (p < zero(p) || p > one(p)) ? NaN : d.location - d.scale*cospi(p)/sinpi(p)
cquantile(d::Cauchy, p::Real) = (p < zero(p) || p > one(p)) ? NaN : d.location + d.scale*cospi(p)/sinpi(p)

mgf(d::Cauchy, t::Real) = NaN
cf(d::Cauchy, t::Real) = exp(im * t * d.location - d.scale * abs(t))

## Sampling
rand(d::Cauchy) = quantile(d,rand())

## Fitting
# Note: this is not a Maximum Likelihood estimator
function fit{T <: Real}(::Type{Cauchy}, x::Array{T})
    l, u = iqr(x)
    Cauchy(median(x), (u - l) / 2.0)
end
