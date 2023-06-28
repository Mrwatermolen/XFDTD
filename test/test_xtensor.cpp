#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xrandom.hpp>

void exampleForDotProduct();
void exampleForCrossProduct();
void exampleForReducerDiff();
void exampleForResize();

int main() {}

void exampleForDotProduct() {
  std::cout << "Xtensor: dot product for two vectors" << std::endl;
  std::vector<size_t> a_shape{3};
  xt::xarray<double> a = xt::random::rand<double>(a_shape, 0.0, 1.0);
  xt::xarray<double> b = xt::random::rand<double>(a_shape, 0.0, 1.0);
  std::cout << a << std::endl;
  std::cout << b << std::endl;
  auto res{0.0};
  std::cout << "a * b: " << xt::linalg::dot(a, b) << std::endl;
  for (auto &&e : a *b) {
    res += e;
  }
  std::cout << "res: " << res << std::endl;
}

void exampleForCrossProduct() {
  std::cout << "Xtensor: cross product for two vectors" << std::endl;
  auto a_shape{std::vector<size_t>{3}};
  auto unit_vec_x{xt::xarray<double>{{1, 0, 0}}};
  auto unit_vec_y{xt::xarray<double>{{0, 1, 0}}};
  std::cout << unit_vec_x << std::endl;
  std::cout << unit_vec_y << std::endl;
  std::cout << "x X y: " << xt::linalg::cross(unit_vec_x, unit_vec_y)
            << std::endl;
}

void exampleForReducerDiff() {
  std::cout << "Xtensor: reducer diff" << std::endl;
  std::cout << "1D" << std::endl;
  auto a{xt::xarray<int>{4, 0, 7, 4, 9, 3, 5, 8, 4, 7, 2, 0}};
  std::cout << a << std::endl;
  auto shifted{xt::view(a, xt::all())};
  std::cout << "Diff: " << xt::diff(shifted) << std::endl;
  std::cout << std::endl;

  std::cout << "2D" << std::endl;
  auto a_shape{std::vector<size_t>{3, 4}};
  a.reshape(a_shape);
  std::cout << a << std::endl;
  // shifted = xt::view(a, xt::all(), xt::all());
  // std::cout << shifted << std::endl;
  std::cout << "Diff: " << xt::diff(a, 1, 0) << std::endl;
}

void exampleForResize() {
  std::cout << "Xtensor: resize" << std::endl;
  xt::xarray<int> a;
  std::cout << a << std::endl;
  std::cout << "shape of a: " << xt::adapt(a.shape()) << std::endl;
  a.resize({1, 4});
  std::cout << a << std::endl;
  std::cout << "shape of a: " << xt::adapt(a.shape()) << std::endl;
  a.resize({4});
  std::cout << a << std::endl;
  std::cout << "shape of a: " << xt::adapt(a.shape()) << std::endl;
}