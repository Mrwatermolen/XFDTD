#include <Eigen/Core>
#include <cstddef>
#include <iostream>
#include <memory>
#include <random>
#include <utility>

#include "shape/cube.h"
#include "shape/shape.h"
#include "shape/sphere.h"

int main() {
  auto cube1{std::make_unique<xfdtd::Cube>(xfdtd::PointVector(0, 0, 0),
                                           xfdtd::PointVector(1, 1, 1))};
  auto sphere1{std::make_unique<xfdtd::Sphere>(xfdtd::PointVector(0, 0, 0), 2)};
  std::shared_ptr<xfdtd::Shape> temp{std::move(cube1->getWrappedBox())};
  auto wrapped_cube1{std::dynamic_pointer_cast<xfdtd::Cube>(temp)};
  temp = std::move(sphere1->getWrappedBox());
  auto wrapped_cube2{std::dynamic_pointer_cast<xfdtd::Cube>(temp)};

  std::cout << "Created cube1: " << static_cast<std::string>(*cube1)
            << std::endl;
  std::cout << "Created sphere1: " << static_cast<std::string>(*sphere1)
            << std::endl;
  std::cout << "Get Box: " << static_cast<std::string>(*wrapped_cube1)
            << std::endl;
  std::cout << "Get Box: " << static_cast<std::string>(*wrapped_cube2)
            << std::endl;

  // Generate random points
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-0.02, 3.0);
  auto dice{[&]() -> double { return distribution(generator); }};
  for (size_t i{0}; i < 100; ++i) {
    auto point{xfdtd::PointVector(dice(), dice(), dice())};
    std::cout << "Point: " << point.transpose() << std::endl;
    std::cout << "isPointInside cube1: " << cube1->isPointInside(point)
              << std::endl;
    std::cout << "isPointInside sphere1: " << sphere1->isPointInside(point)
              << std::endl;
  }
}
