#include "LaplaceDistribution.hpp"
#include <cmath>

namespace ptm {

LaplaceDistribution::LaplaceDistribution(double mu, double b) : mu_(mu), b_(b) {
}

double LaplaceDistribution::Pdf(double x) const {
  return (1.0 / (2.0 * b_)) * std::exp(-std::abs(x - mu_) / b_);
}

double LaplaceDistribution::Cdf(double x) const {
  if (x < mu_)
    return 0.5 * std::exp((x - mu_) / b_);
  return 1.0 - 0.5 * std::exp(-(x - mu_) / b_);
}

double LaplaceDistribution::Sample(std::mt19937& rng) const {
  std::uniform_real_distribution<double> dist(-0.5, 0.5);
  double u = dist(rng);
  return mu_ - b_ * (u < 0 ? -1.0 : 1.0) * std::log(1.0 - 2.0 * std::abs(u));
}

double LaplaceDistribution::TheoreticalMean() const {
  return mu_;
}
double LaplaceDistribution::TheoreticalVariance() const {
  return 2.0 * b_ * b_;
}

} // namespace ptm
