#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#include "waveform/gaussian_waveform.h"

int main() {
  constexpr int nc{20};
  constexpr double c{3.0e8};
  constexpr double dx{1e-3};
  constexpr double dt{dx / (c * 2)};
  constexpr double tau{nc * dx / (2 * c)};
  constexpr double t_0{4.5 * tau};
  auto gaussian_waveform{xfdtd::GaussianWaveform(1, tau, t_0)};
  size_t time_steps{1000};
  auto time_array{std::vector<double>{}};
  time_array.resize(time_steps);
  for (size_t i{0}; i < time_steps; ++i) {
    time_array[i] = dt * i;
  }
  gaussian_waveform.init(time_array);
  std::ofstream data_file("test_source.dat");
  if (!data_file.is_open()) {
    std::cout << "Error opening file" << std::endl;
    return 1;
  }

  for (size_t i{0}; i < time_steps; ++i) {
    data_file << time_array[i] << " " << gaussian_waveform.getValue(i)
              << std::endl;
  }
}
