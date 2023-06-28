#include <fftw3.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <tuple>
#include <vector>
#include <xtensor-fftw/basic.hpp>   // rfft, irfft
#include <xtensor-fftw/helper.hpp>  // rfftscale
#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xbuilder.hpp>  // xt::arange
#include <xtensor/xio.hpp>
#include <xtensor/xmanipulation.hpp>
#include <xtensor/xmath.hpp>  // xt::sin, cos
#include <xtensor/xview.hpp>

#include "util/constant.h"
#include "xtensor-fftw/basic_double.hpp"

#define NUM_POINTS 64
#define REAL 0
#define IMAG 1

template <typename T>
xt::xarray<std::complex<T>> fftshift(const xt::xarray<std::complex<T>>& input) {
  auto shape = input.shape();
  std::vector<size_t> shift(shape.size());
  for (size_t i = 0; i < shape.size(); ++i) {
    shift[i] = shape[i] / 2;
  }
  return xt::roll(input, shift);
}

int main() {
  // generate a sinusoid field
  double dx = M_PI / 100;
  xt::xarray<double> x = xt::arange(0., 2 * M_PI, dx);
  xt::xarray<double> sin = xt::sin(x);

  // transform to Fourier space
  auto sin_fs = xt::fftw::rfft(sin);

  // multiply by i*k
  std::complex<double> i{0, 1};
  auto k = xt::fftw::rfftscale<double>(sin.shape()[0], dx);
  xt::xarray<std::complex<double>> sin_derivative_fs = xt::eval(i * k * sin_fs);

  // transform back to normal space
  auto sin_derivative = xt::fftw::irfft(sin_derivative_fs);

  // std::cout << "x:              " << x << std::endl;
  // std::cout << "sin:            " << sin << std::endl;
  // std::cout << "cos:            " << xt::cos(x) << std::endl;
  // std::cout << "sin_derivative: " << sin_derivative << std::endl;

  std::vector<size_t> a_shape{2, 3, 4};
  auto a{xt::xarray<double>(a_shape)};
  double c{0};
  for (auto&& e : a) {
    e = c++;
  }
  std::cout << a << std::endl;
  std::cout << xt::xarray<double>(xt::view(a, xt::all(), 2, -1)) << std::endl;

  // Create a 3D array with shape (3, 4, 4)
  xt::xarray<int> data = xt::zeros<int>({3, 4, 4});

  std::cout << "Original 3D array:" << std::endl;
  std::cout << data << std::endl;

  // Modify the second 2D slice along the first axis (axis 0)
  xt::view(data, 1, xt::all(), xt::all()) = xt::arange<int>(16).reshape({4, 4});
  auto new_data = xt::ones<int>({4, 4});
  xt::view(data, 2, xt::all(), xt::all()) = new_data;

  std::cout << "Modified 3D array:" << std::endl;
  std::cout << data << std::endl;

  std::cout << "xt::concatenate" << std::endl;
  a = xt::ones<int>(std::vector<size_t>{3, 3, 3});
  auto b = xt::zeros<int>(std::vector<size_t>{3, 3, 2});
  std::cout << a << std::endl;
  std::cout << b << std::endl;

  std::cout << xt::concatenate(xt::xtuple(a, b), 2) << std::endl;

  std::cout << xt::view(data, xt::all(), 0, 0).shape(0)
            << xt::view(data, xt::all(), 0, 0).shape(1) << std::endl;

  return 0;
}