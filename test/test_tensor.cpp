#include <Eigen/Core>
#include <array>
#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>

int main() {
  using Eigen::Tensor;
  Tensor<double, 3> ex(2, 2, 2);
  Tensor<double, 3> ey(2, 2, 2);
  ex.setZero();
  ey.setZero();

  std::cout << ex << std::endl;
  std::cout << ey << std::endl;

  ex(0, 0, 0) = 1;
  ex(0, 0, 1) = 2;
  ex(0, 1, 0) = 3;
  ex(0, 1, 1) = 4;
  ex(1, 0, 0) = 5;
  ex(1, 0, 1) = 6;
  ex(1, 1, 0) = 7;
  ex(1, 1, 1) = 8;
  ey(0, 0, 0) = 1;
  ey(0, 0, 1) = 2;
  ey(0, 1, 0) = 3;
  ey(0, 1, 1) = 4;
  ey(1, 0, 0) = 5;
  ey(1, 0, 1) = 6;
  ey(1, 1, 0) = 7;
  ey(1, 1, 1) = 8;

  std::cout << ex << std::endl;
  std::cout << ey << std::endl;

  std::cout << ex * ey << std::endl;

  auto res{ex * ey};
  auto r{ex(1, 1, 1) * ey(0, 1, 1)};
  std::cout << res << std::endl;
  std::cout << r << std::endl;

  ex(1, 1, 1) = ex(1, 1, 1) * ey(0, 1, 1);
  std::cout << ex(1, 1, 1) << std::endl;

  Tensor<double, 3> tensor(3, 4, 5);
  tensor.setRandom();

  Eigen::array<Eigen::Index, 3> dims {tensor.dimensions()};
  std::array<int, 3> offset = {0,0,0};
  std::array<int, 3> extent = {1,4,5};
  std::array<int, 3> shape = {3,4};
  
  tensor.slice(offset, extent).reshape(shape);
  return 0;
}
