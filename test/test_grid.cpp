#include <grid/grid_space.h>
#include <shape/cube.h>

#include <cassert>
#include <cstddef>
#include <iostream>
#include <random>

#include "util/float_compare.h"
#include "util/type_define.h"

void testGridSpace1D(int times) {
  using xfdtd::Cube;
  using xfdtd::Grid;
  using xfdtd::GridSpace;
  using xfdtd::PointVector;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<std::size_t> dis(0, -1);
  auto random_seed{dis(gen)};
  std::clog << "Test GridSpace 1D with random seed: " << random_seed << "\n";

  std::default_random_engine e(random_seed);
  std::uniform_real_distribution<double> u(0, 1);
  for (int i{0}; i < times; ++i) {
    std::clog << "Test GridSpace 1D Number: " << i << "\n";
    auto origin_point{PointVector{u(e), u(e), u(e)}};
    auto size_x{u(e)};
    auto size{PointVector{size_x, size_x, size_x}};
    auto cube{std::make_unique<Cube>(origin_point, size)};
    auto dl{cube->getSize()[2] / 20.0};
    auto grid_space{std::make_shared<GridSpace>(dl, dl, dl)};
    grid_space->calculateSpaceSize(cube->getWrappedBox().get());
    grid_space->generateUniformGridSpace();
    grid_space->tellMeOk();

    assert(grid_space->getDimension() == GridSpace::Dimension::THREE_DIMENSION);

    auto nx{20};
    auto ny{20};
    auto nz{20};

    assert(grid_space->getGridNumX() == nx);
    assert(grid_space->getGridNumY() == ny);
    assert(grid_space->getGridNumZ() == nz);

    auto dx{dl};
    auto dy{dl};
    auto dz{dl};

    assert(grid_space->getGridBaseSizeX() == dx);
    assert(grid_space->getGridBaseSizeY() == dy);
    assert(grid_space->getGridBaseSizeZ() == dz);

    assert(grid_space->getGridSizeMinX() == dz);
    assert(grid_space->getGridSizeMinY() == dy);
    assert(grid_space->getGridSizeMinZ() == dx);

    auto& grid_size_array_x{grid_space->getGridSizeArrayX()};
    auto& grid_size_array_y{grid_space->getGridSizeArrayY()};
    auto& grid_size_array_z{grid_space->getGridSizeArrayZ()};
    for (int k{0}; k < nz; ++k) {
      assert(xfdtd::isEqual(grid_size_array_x(k), dx));
      assert(xfdtd::isEqual(grid_size_array_y(k), dy));
      assert(xfdtd::isEqual(dz, grid_size_array_z(k)));
    }

    auto& grid_size_array_ex{grid_space->getGridSizeArrayEX()};

    for (int k{0}; k < nz; ++k) {
      assert(xfdtd::isEqual(grid_size_array_ex(k), dx));
    }

    auto& grid_size_array_ey{grid_space->getGridSizeArrayEY()};

    auto& grid_size_array_ez{grid_space->getGridSizeArrayEZ()};

    for (int k{0}; k < nz; ++k) {
      assert(xfdtd::isEqual(grid_size_array_ex(k), dx));
      assert(xfdtd::isEqual(grid_size_array_ey(k), dy));
      assert(xfdtd::isEqual(grid_size_array_ez(k), dz));
    }

    for(auto k{0}; k < nz; ++k) {
      assert(xfdtd::isEqual(grid_space->getGridSizeHX(k), grid_size_array_ex[k]));
      assert(xfdtd::isEqual(grid_space->getGridSizeHY(k), grid_size_array_ey[k]));
      assert(xfdtd::isEqual(grid_space->getGridSizeHZ(k), grid_size_array_ez[k]));
    }
  }
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cout << "argc: " << argc << "\n";
    std::cerr << "Usage: " << argv[0] << " <test name> <test times>\n";
    return 1;
  }

  auto test_name{std::string(argv[1])};
  auto test_times{std::stoi(argv[2])};

  if (test_name == "GridSpace1D") {
    testGridSpace1D(test_times);
  } else {
    std::cerr << "Unknown test name: " << test_name << "\n";
    return 1;
  }
}