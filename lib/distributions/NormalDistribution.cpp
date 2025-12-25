#include "NormalDistribution.hpp"
#include <cmath>

namespace ptm {

NormalDistribution::NormalDistribution(double mean, double stddev) : mean_(mean), stddev_(stddev) {
}

double NormalDistribution::Pdf(double x) const {
  static const double inv_sqrt_2pi = 1.0 / std::sqrt(2.0 * M_PI);
  double z = (x - mean_) / stddev_;
  return (inv_sqrt_2pi / stddev_) * std::exp(-0.5 * z * z);
}

double NormalDistribution::Cdf(double x) const {
  return 0.5 * (1.0 + std::erf((x - mean_) / (stddev_ * std::numbers::sqrt2)));
}

double NormalDistribution::Sample(std::mt19937& rng) const {
  std::normal_distribution<double> dist(mean_, stddev_);
  return dist(rng);
}

double NormalDistribution::TheoreticalMean() const {
  return mean_;
}
double NormalDistribution::TheoreticalVariance() const {
  return stddev_ * stddev_;
}
double NormalDistribution::GetMean() const {
  return mean_;
}
double NormalDistribution::GetStddev() const {
  return stddev_;
}

} // namespace ptm
