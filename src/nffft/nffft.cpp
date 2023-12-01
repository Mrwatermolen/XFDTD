#include "nffft/nffft.h"

#include <string>
#include <utility>
#include <xtensor-fftw/basic_double.hpp>
#include <xtensor-fftw/helper.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>

#include "util/type_define.h"

namespace xfdtd {

NFFFT::NFFFT(SpatialIndex distance_x, SpatialIndex distance_y,
             SpatialIndex distance_z, std::string output_dir_path)
    : _distance_x{distance_x},
      _distance_y{distance_y},
      _distance_z{distance_z},
      _output_dir_path{std::move(output_dir_path)} {}

void NFFFT::defaultInit(std::shared_ptr<const GridSpace> grid_space,
                        std::shared_ptr<const FDTDBasicCoff> fdtd_basic_coff,
                        std::shared_ptr<const EMF> emf) {
  _grid_space = std::move(grid_space);
  _fdtd_basic_coff = std::move(fdtd_basic_coff);
  _emf = std::move(emf);

  _output_box =
      std::make_unique<GridBox>(_distance_x, _distance_y, _distance_z,
                                _grid_space->getGridNumX() - 2 * _distance_x,
                                _grid_space->getGridNumY() - 2 * _distance_y,
                                _grid_space->getGridNumZ() - 2 * _distance_z);
}

}  // namespace xfdtd
