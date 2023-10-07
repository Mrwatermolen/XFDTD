#include "nffft/nffft.h"

#include <utility>
#include <xtensor-fftw/basic_double.hpp>
#include <xtensor-fftw/helper.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>

#include "util/type_define.h"

namespace xfdtd {

NFFFT::NFFFT(SpatialIndex distance_x, SpatialIndex distance_y,
             SpatialIndex distance_z, std::filesystem::path output_dir_path)
    : _distance_x{distance_x},
      _distance_y{distance_y},
      _distance_z{distance_z},
      _output_dir_path{std::move(output_dir_path)} {}

void NFFFT::defaultInit(std::unique_ptr<GridBox> output_box,
                        std::shared_ptr<EMF> emf, size_t total_time_steps,
                        double dt, double dx, double dy, double dz) {
  if (output_box == nullptr) {
    throw std::runtime_error("Output box instance is not set.");
  }
  if (emf == nullptr) {
    throw std::runtime_error("EMF instance is not set.");
  }

  _output_box = std::move(output_box);
  _emf = std::move(emf);
  _total_time_steps = total_time_steps;
  _dt = dt;
  _dx = dx;
  _dy = dy;
  _dz = dz;
}

}  // namespace xfdtd
