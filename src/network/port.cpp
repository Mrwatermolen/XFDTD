#include "network/port.h"

#include <complex>
#include <utility>

#include "util/dft.h"

namespace xfdtd {

Port::Port(int index, std::complex<double> impedance, bool is_activated_source,
           std::shared_ptr<xfdtd::VoltageMonitor> voltage,
           std::shared_ptr<xfdtd::CurrentMonitor> current)
    : _index{index},
      _impedance{impedance},
      _is_activated_source{is_activated_source},
      _voltage{std::move(voltage)},
      _current{std::move(current)} {}

void Port::init(const std::shared_ptr<const FDTDBasicCoff> &fdtd_basic_coff,
                const std::shared_ptr<const GridSpace> &grid_space,
                const std::shared_ptr<const EMF> &emf) {
  _dt = fdtd_basic_coff->getDt();
}

void Port::calculateSParameter(const xt::xarray<double> &_frequencies) {
  auto dt{_dt};
  auto z{_impedance};
  auto k{std::sqrt(std::real(z))};

  auto v{dft(_voltage->getValueArray(), dt, _frequencies)};
  auto i{dft(_current->getValueArray(), dt, _frequencies)};
  _a = 0.5 * (v + z * i) / k;
  _b = 0.5 * (v - std::conj(z) * i) / k;
}

int Port::getPortIndex() const { return _index; }

bool Port::isActivatedSource() const { return _is_activated_source; }

xt::xarray<std::complex<double>> Port::getA() const { return _a; }

xt::xarray<std::complex<double>> Port::getB() const { return _b; }

}  // namespace xfdtd
