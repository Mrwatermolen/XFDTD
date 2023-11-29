#include "monitor/current_monitor.h"

#include <filesystem>
#include <xtensor/xnpy.hpp>

#include "util/type_define.h"

namespace xfdtd {

CurrentMonitor::CurrentMonitor(std::unique_ptr<Cube> shape,
                               Orientation orientation,
                               std::filesystem::path output_dir_path,
                               std::string output_file_name)
    : Monitor{std::move(shape), std::move(output_dir_path),
              std::move(output_file_name)},
      _orientation{orientation} {}

std::unique_ptr<Monitor> CurrentMonitor::clone() const {
  return std::make_unique<CurrentMonitor>(*this);
}

void CurrentMonitor::init(
    const std::shared_ptr<const FDTDBasicCoff> &fdtd_basic_coff,
    const std::shared_ptr<const GridSpace> &grid_space,
    const std::shared_ptr<const EMF> &emf) {
  defaultInit(fdtd_basic_coff, grid_space, emf);
  auto grid_box{grid_space->getGridBox(getShape())};

  _is = grid_box->getGridStartIndexX();
  _ie = grid_box->getGridEndIndexX();
  _js = grid_box->getGridStartIndexY();
  _je = grid_box->getGridEndIndexY();
  _ks = grid_box->getGridStartIndexZ();
  _ke = grid_box->getGridEndIndexZ();

  if (_orientation == Orientation::XN or _orientation == Orientation::XP) {
    _da = xt::view(grid_space->getGridSizeArrayEY(), xt::range(_js, _je + 1));
    _db = xt::view(grid_space->getGridSizeArrayEZ(), xt::range(_ks, _ke + 1));
    if (_orientation == Orientation::XN) {
      _positive = -1;
    } else {
      _positive = 1;
    }
  }

  if (_orientation == Orientation::YN or _orientation == Orientation::YP) {
    _da = xt::view(grid_space->getGridSizeArrayEZ(), xt::range(_ks, _ke + 1));
    _db = xt::view(grid_space->getGridSizeArrayEX(), xt::range(_is, _ie + 1));
    if (_orientation == Orientation::YN) {
      _positive = -1;
    } else {
      _positive = 1;
    }
  }

  if (_orientation == Orientation::ZN or _orientation == Orientation::ZP) {
    _da = xt::view(grid_space->getGridSizeArrayEX(), xt::range(_is, _ie + 1));
    _db = xt::view(grid_space->getGridSizeArrayEY(), xt::range(_js, _je + 1));
    if (_orientation == Orientation::ZN) {
      _positive = -1;
    } else {
      _positive = 1;
    }
  }

  _value = xt::zeros<double>({getFDTDBasicCoffInstance()->getTotalTimeStep()});
}

void CurrentMonitor::update() {
  auto emf{getEMFInstance()};
  auto t{getFDTDBasicCoffInstance()->getCurrentTimeStep()};
  double current{0};
  double integral_x{0};
  double integral_y{0};
  double integral_z{0};

  if (_orientation == Orientation::XN or _orientation == Orientation::XP) {
    for (size_t j{_js}; j <= _je; ++j) {
      integral_y += (emf->getHy(_is, j, _ks - 1) - emf->getHy(_is, j, _ke)) *
                    _da(j - _js);
    }
    for (size_t k{_ks}; k <= _ke; ++k) {
      integral_z += (emf->getHz(_is, _je, k) - emf->getHz(_is, _js - 1, k)) *
                    _db(k - _ks);
    }
    _value(t) = _positive * (integral_z + integral_y);
    return;
  }

  if (_orientation == Orientation::YN || _orientation == Orientation::YP) {
    for (size_t k{_ks}; k <= _ke; ++k) {
      integral_z += (emf->getHz(_is - 1, _js, k) - emf->getHz(_ie, _js, k)) *
                    _da(k - _ks);
    }
    for (size_t i{_is}; i <= _ie; ++i) {
      integral_x += (emf->getHx(i, _js, _ke) - emf->getHx(i, _js, _ks - 1)) *
                    _db(i - _is);
    }
    _value(t) = _positive * (integral_z + integral_x);
    return;
  }
  if (_orientation == Orientation::ZN || _orientation == Orientation::ZP) {
    for (size_t i{_is}; i <= _ie; ++i) {
      integral_x += (emf->getHx(i, _js - 1, _ks) - emf->getHx(i, _je, _ks));
    }
    for (size_t j{_js}; j <= _je; ++j) {
      integral_y += (emf->getHy(_ie, j, _ks) - emf->getHy(_is - 1, j, _ks));
    }
    _value(t) = _positive * (integral_x + integral_y);
    return;
  }
}

void CurrentMonitor::outputData() {
  if (!std::filesystem::exists(getOutputPath()) &&
      !std::filesystem::is_directory(getOutputPath())) {
    try {
      std::filesystem::create_directories(getOutputPath());
    } catch (std::exception e) {
      std::cerr << "Error: cannot create directory " << getOutputPath() << '\n';
      return;
    }
  }

  auto time_array{xt::linspace<double>(
      0.5 * getFDTDBasicCoffInstance()->getDt(),
      (getFDTDBasicCoffInstance()->getTotalTimeStep() + 0.5) *
          getFDTDBasicCoffInstance()->getDt(),
      getFDTDBasicCoffInstance()->getTotalTimeStep())};
  auto temp{xt::stack(xt::xtuple(time_array, _value))};

  xt::dump_npy(getOutputPath() / getOutputFileName(), temp);
}

double CurrentMonitor::getValue(size_t i) const { return _value(i); }

}  // namespace xfdtd
