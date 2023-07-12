#ifndef _DFT_H_
#define _DFT_H_

#include <complex>
#include <cstddef>
#include <xtensor/xarray.hpp>

#include "util/constant.h"

namespace xfdtd {
inline xt::xarray<std::complex<double>> dft(
    const xt::xarray<double> &time_domain_data, double dt,
    const xt::xarray<double> &frequencies) {
  using namespace std::complex_literals;
  xt::xarray<std::complex<double>> res;
  res.resize({frequencies.size()});
  for (size_t i{0}; i < frequencies.size(); ++i) {
    std::complex<double> sum{0.0, 0.0};
    for (size_t j{0}; j < time_domain_data.size(); ++j) {
      auto t{j * dt};
      sum += time_domain_data(j) *
             std::exp(-2.0 * 1i * constant::PI * frequencies(i) * t);
    }
    res(i) = sum * dt;
  }
  return res;
}
}  // namespace xfdtd

#endif  // _DFT_H_