#include <memory>
#include <utility>
#include <vector>

#include "boundary/perfect_match_layer.h"
#include "lumped_element/lumped_element.h"
#include "lumped_element/voltage_source.h"
#include "material/material.h"
#include "monitor/current_monitor.h"
#include "monitor/voltage_monitor.h"
#include "network/port.h"
#include "nffft/nffft_fd.h"
#include "object/object.h"
#include "object/object_plane.h"
#include "shape/cube.h"
#include "simulation/simulation.h"
#include "util/constant.h"
#include "util/type_define.h"
#include "waveform/gaussian_waveform.h"

// TODO(franzero): can't run this test case. some bugs in the code.
void testInvertedFAntenna() {
  using namespace xfdtd;

  constexpr double cell_size{2e-4};

  Cube cube{{0, 0, 0}, {1, 1, 1}};
  // define material
  auto air{Material{"air", 1, 1, 0, 0}};
  auto pec{Material{"pec", 1, 1, 1e10, 0}};
  auto substrate{Material{"substrate", 2.2, 1, 0, 0}};

  // define object
  // domain
  auto domain_origin_point{PointVector{0, -5 * cell_size, -5 * cell_size}};
  auto domain_size{PointVector{6 * cell_size, 40 * cell_size, 40 * cell_size}};
  auto domain_shape{Cube{domain_origin_point, domain_size}};
  auto domain_object{Object{"domain", std::make_unique<Cube>(domain_shape),
                            std::make_unique<Material>(air)}};
  // substrate
  auto substrate_origin_point{PointVector{-0.787e-3, 0, 0}};
  auto substrate_size{PointVector{0.787e-3, 40e-3, 40e-3}};
  auto substrate_shape{Cube{substrate_origin_point, substrate_size}};
  auto substrate_object{Object{"substrate",
                               std::make_unique<Cube>(substrate_shape),
                               std::make_unique<Material>(substrate)}};

  //
  auto cube_0_origin_point{PointVector{0, 0, 24e-3}};
  auto cube_0_size{PointVector{0, 28.4e-3, 2.4e-3}};
  auto cube_0_shape{Cube{cube_0_origin_point, cube_0_size}};
  auto cube_0_object{ObjectPlane{"cube_0", std::make_unique<Cube>(cube_0_shape),
                                 std::make_unique<Material>(pec)}};

  auto cube_1_origin_point{PointVector{0, 16e-3, 30e-3}};
  auto cube_1_size{PointVector{0, 12.4e-3, 2.4e-3}};
  auto cube_1_shape{Cube{cube_1_origin_point, cube_1_size}};
  auto cube_1_object{ObjectPlane{"cube_1", std::make_unique<Cube>(cube_1_shape),
                                 std::make_unique<Material>(pec)}};

  auto cube_2_origin_point{PointVector{0, 26e-3, 8.4e-3}};
  auto cube_2_size{PointVector{0, 2.4e-3, 24.4e-3}};
  auto cube_2_shape{Cube{cube_2_origin_point, cube_2_size}};
  auto cube_2_object{ObjectPlane{"cube_2", std::make_unique<Cube>(cube_2_shape),
                                 std::make_unique<Material>(pec)}};

  auto cube_3_origin_point{PointVector{0, 20.8e-3, 16e-3}};
  auto cube_3_size{PointVector{0, 2.4e-3, 16.4e-3}};
  auto cube_3_shape{Cube{cube_3_origin_point, cube_3_size}};
  auto cube_3_object{ObjectPlane{"cube_3", std::make_unique<Cube>(cube_3_shape),
                                 std::make_unique<Material>(pec)}};

  auto cube_4_origin_point{PointVector{-0.787e-3, 15.6e-3, 24e-3}};
  auto cube_4_size{PointVector{0.787e-3, 0, 2.4e-3}};
  auto cube_4_shape{Cube{cube_4_origin_point, cube_4_size}};
  auto cube_4_object{ObjectPlane{"cube_4", std::make_unique<Cube>(cube_4_shape),
                                 std::make_unique<Material>(pec)}};

  auto cube_5_origin_point{PointVector{-0.787e-3, 0, 0}};
  auto cube_5_size{PointVector{0, 16e-3, 40e-3}};
  auto cube_5_shape{Cube{cube_5_origin_point, cube_5_size}};
  auto cube_5_object{ObjectPlane{"cube_5", std::make_unique<Cube>(cube_5_shape),
                                 std::make_unique<Material>(pec)}};
  std::vector<std::shared_ptr<Object>>objects;
  objects.emplace_back(std::make_unique<Object>(domain_object));
  objects.emplace_back(std::make_unique<Object>(substrate_object));
  objects.emplace_back(std::make_unique<ObjectPlane>(cube_0_object));
  objects.emplace_back(std::make_unique<ObjectPlane>(cube_1_object));
  objects.emplace_back(std::make_unique<ObjectPlane>(cube_2_object));
  objects.emplace_back(std::make_unique<ObjectPlane>(cube_3_object));
  objects.emplace_back(std::make_unique<ObjectPlane>(cube_4_object));
  objects.emplace_back(std::make_unique<ObjectPlane>(cube_5_object));

  // lumped element
  auto voltage_source_origin_point{PointVector{-0.787e-3, 0, 24e-3}};
  auto voltage_source_size{PointVector{0, 1 * cell_size, 2.4e-3}};
  auto voltage_source_shape{
      Cube{voltage_source_origin_point, voltage_source_size}};
  constexpr double l_min{20 * cell_size};
  constexpr double tau{l_min / 6e8};
  constexpr double t_0{4.5 * tau};
  auto waveform{GaussianWaveform{1, tau, t_0}};
  auto voltage_source{VoltageSource{
      std::make_unique<Cube>(voltage_source_shape), Orientation ::XP, 50,
      std::make_unique<GaussianWaveform>(waveform)}};

  std::vector<std::shared_ptr<LumpedElement>> lumped_elements;
  lumped_elements.emplace_back(std::make_shared<VoltageSource>(voltage_source));

  // monitor
  auto v1{std::make_shared<VoltageMonitor>(
      std::make_unique<Cube>(PointVector{-0.787e-3, 0, 24e-3},
                             PointVector{0, 1 * cell_size, 2.4e-3}),
      Orientation::XP, "./visualizing_data/data/inverted_f_antenna", "v1.npy")};
  auto c1{std::make_shared<CurrentMonitor>(
      std::make_unique<Cube>(PointVector{-0.39e-3, 0, 24e-33},
                             PointVector{1 * cell_size, 1 * cell_size, 2.4e-3}),
      Orientation::XP, "./visualizing_data/data/inverted_f_antenna", "c1.npy")};

  std::vector<std::shared_ptr<Monitor>> monitors;
  monitors.emplace_back(v1);
  monitors.emplace_back(c1);

  // network
  auto port_1{Port{1, 50, true, v1, c1}};
  std::vector<std::unique_ptr<Port>> ports;
  ports.emplace_back(std::make_unique<Port>(std::move(port_1)));
  auto network{
      std::make_unique<Network>(std::move(ports), xt::arange(2e7, 1e10, 2e7),
                                "./visualizing_data/data/inverted_f_antenna")};

  // boundary
  std::vector<std::shared_ptr<Boundary>> boundaries;
  boundaries.emplace_back(std::make_shared<PML>(Orientation::XN, 8));
  boundaries.emplace_back(std::make_shared<PML>(Orientation::XP, 8));
  boundaries.emplace_back(std::make_shared<PML>(Orientation::YN, 8));
  boundaries.emplace_back(std::make_shared<PML>(Orientation::YP, 8));
  boundaries.emplace_back(std::make_shared<PML>(Orientation::ZN, 8));
  boundaries.emplace_back(std::make_shared<PML>(Orientation::ZP, 8));

  // FDTD
  Simulation simulation{
      cell_size,
      std::move(objects),
      std::move(lumped_elements),
      std::move(boundaries),
      std::move(monitors),
      std::move(network),
  };
  simulation.addNFFFT(std::make_unique<NffftFd>(
      13, 13, 13, xt::xarray<double>({2.4e9, 5.8e9}), xt::xarray<double>({0}),
      xt::linspace(0.0, constant::PI * 2, 360),
      "./visualizing_data/data/inverted_f_antenna"));

  simulation.run(2000);
}

int main() {
  std::cerr << "Error: This case (InvertedFAntenna) can't run currently.\n";
  return 0;
}