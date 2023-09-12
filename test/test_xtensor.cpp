#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor.hpp>
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
void exampleForStack();
void exampleForMeshgrid();
void exampleForFind();

int main() { exampleForFind(); }

void exampleForDotProduct() {
  std::cout << "Xtensor: dot product for two vectors" << std::endl;
  std::vector<size_t> a_shape{3};
  xt::xarray<double> a = xt::random::rand<double>(a_shape, 0.0, 1.0);
  xt::xarray<double> b = xt::random::rand<double>(a_shape, 0.0, 1.0);
  std::cout << a << std::endl;
  std::cout << b << std::endl;
  auto res{0.0};
  std::cout << "a * b: " << xt::linalg::dot(a, b) << std::endl;
  for (auto&& e : a* b) {
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

void exampleForStack() {
  std::cout << "Xtensor: stack" << std::endl;
  size_t l{0};
  size_t r{4};
  auto i_range{xt::arange<double>(l, r)};
  auto j_range{xt::arange<double>(l, r)};
  auto k_range{xt::arange<double>(l, r)};
  j_range + 0.5;
  auto unit_vec{xt::xtensor_fixed<double, xt::xshape<3>>{
      sin(M_PI / 4) * cos(M_PI / 4), sin(M_PI / 4) * sin(M_PI / 4),
      cos(M_PI / 4)}};
  auto matrix{xt::stack(xt::xtuple(i_range, j_range + 0.5, k_range + 1), 1)};
  std::cout << matrix << std::endl;
  std::cout << "Vector:" << unit_vec << std::endl;
  std::cout << xt::linalg::dot(matrix, unit_vec);
}

void exampleForMeshgrid() {
  auto [x_matrix, y_matrix, z_matrix] =
      xt::meshgrid(xt::linspace(0.0, 1.0, 2), xt::arange<double>(1, 1, 1),
                   xt::linspace(0.0, 1.0, 2));

  std::cout << "X:\n" << x_matrix << std::endl;
  std::cout << "Y:\n" << y_matrix << std::endl;
  std::cout << "Z:\n" << z_matrix << std::endl;

  auto r{xt::xarray<double>{0.5, 0.5, 0.7}};
  std::cout << "R:\n" << r << std::endl;

  auto x_r{x_matrix * r(0)};
  auto y_r{y_matrix * r(1)};
  auto z_r{z_matrix * r(2)};

  std::cout << "X_r:\n" << x_r << std::endl;
  std::cout << "Y_r:\n" << y_r << std::endl;
  std::cout << "Z_r:\n" << z_r << std::endl;

  xt::xarray<double> sum{x_r + y_r + z_r};
  std::cout << "Sum:\n" << xt::transpose(sum) << std::endl;

  auto matrix =
      xt::stack(xt::xtuple(x_matrix.reshape({-1}), y_matrix.reshape({-1}),
                           z_matrix.reshape({-1})),
                1);

  std::cout << "Matrix:\n" << matrix << std::endl;

  auto res{xt::linalg::dot(matrix, r)};
  std::cout << "Result:\n" << res << std::endl;
}

void exampleForFind() {
  // 创建一个 3x3 的 xtensor 数组
  xt::xarray<int> arr{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  xt::xarray<int> coff{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  // 使用一个 n 维的 std::vector 作为索引来修改元素
  std::vector<std::size_t> index = {1, 2};
  xt::xarray<bool> mask = arr > 3;
  xt::filter(coff, mask) = 0;

  // 打印修改后的数组
  std::cout << arr << std::endl;
  std::cout << coff << std::endl;
  for (auto e : arr) {
    std::cout << e << "\t";
  }
  std::cout << std::endl;
}