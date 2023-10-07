#include "monitor/voltage_monitor.h"

#include <utility>
#include <xtensor/xindex_view.hpp>
#include <xtensor/xnpy.hpp>

#include "util/type_define.h"

namespace xfdtd {

VoltageMonitor::VoltageMonitor(std::unique_ptr<Cube> shape,
                               Orientation orientation,
                               std::filesystem::path output_dir_path,
                               std::string output_file_name)
    : Monitor{std::move(shape), std::move(output_dir_path),
              std::move(output_file_name)},
      _orientation{orientation} {}

std::unique_ptr<Monitor> VoltageMonitor::clone() const {
  return std::make_unique<VoltageMonitor>(*this);
}

void VoltageMonitor::init(const std::shared_ptr<const FDTDBasicCoff> &fdtd_basic_coff,
                const std::shared_ptr<const GridSpace> &grid_space,
                const std::shared_ptr<const EMF> &emf) {
  defaultInit(fdtd_basic_coff, grid_space, emf);

  auto grid_box{grid_space->getGridBox(getShape())};
  _is = grid_box->getGridOriginIndexX();
  _js = grid_box->getGridOriginIndexY();
  _ks = grid_box->getGridOriginIndexZ();
  _ie = grid_box->getGridEndIndexX();
  _je = grid_box->getGridEndIndexY();
  _ke = grid_box->getGridEndIndexZ();

  if (_orientation == Orientation::XN || _orientation == Orientation::XP) {
    _dc = fdtd_basic_coff->getDx();
    _coff = -_dc / ((_je - _js) * (_ke - _ks));
  }
  if (_orientation == Orientation::YN || _orientation == Orientation::YP) {
    _dc = fdtd_basic_coff->getDy();
    _coff = -_dc / ((_ie - _is) * (_ke - _ks));
  }
  if (_orientation == Orientation::ZN || _orientation == Orientation::ZP) {
    _dc = fdtd_basic_coff->getDz();
    _coff = -_dc / ((_ie - _is) * (_je - _js));
    if (_orientation == Orientation::ZN) {
      _coff *= -1;
    }
  }

  _value.resize({getFDTDBasicCoffInstance()->getTotalTimeStep()});
}

void VoltageMonitor::update() {
  auto current_time_step{getFDTDBasicCoffInstance()->getCurrentTimeStep()};
  const auto emf{getEMFInstance()};

  if (_orientation == Orientation::XN || _orientation == Orientation::XP) {
    for (size_t i{_is}; i < _ie; ++i) {
      for (size_t j{_js}; j < _je; ++j) {
        for (size_t k{_ks}; k < _ke; ++k) {
          _value(current_time_step) += _coff * (emf->getEx(i, j, k));
        }
      }
    }
  }

  if (_orientation == Orientation::YN || _orientation == Orientation::YP) {
    for (size_t i{_is}; i < _ie; ++i) {
      for (size_t j{_js}; j < _je; ++j) {
        for (size_t k{_ks}; k < _ke; ++k) {
          _value(current_time_step) += _coff * (emf->getEy(i, j, k));
        }
      }
    }
  }

  if (_orientation == Orientation::ZN || _orientation == Orientation::ZP) {
    for (size_t i{_is}; i < _ie; ++i) {
      for (size_t j{_js}; j < _je; ++j) {
        for (size_t k{_ks}; k < _ke; ++k) {
          _value(current_time_step) += _coff * (emf->getEz(i, j, k));
        }
      }
    }
  }
}

void VoltageMonitor::outputData() {
  if (!std::filesystem::exists(getOutputPath()) &&
      !std::filesystem::is_directory(getOutputPath())) {
    try {
      std::filesystem::create_directory(getOutputPath());
    } catch (std::exception e) {
      std::cerr << "Error: cannot create directory " << getOutputPath() << '\n';
      return;
    }
  }

  auto time_array{
      xt::linspace<double>(0.0,
                           getFDTDBasicCoffInstance()->getTotalTimeStep() *
                               getFDTDBasicCoffInstance()->getDt(),
                           getFDTDBasicCoffInstance()->getTotalTimeStep())};
  auto temp{xt::stack(xt::xtuple(time_array, _value))};

  xt::dump_npy(getOutputPath() / getOutputFileName(), temp);
}

double VoltageMonitor::getValue(size_t i) const { return _value(i); }

}  // namespace xfdtd
