#ifndef _XFDTD_PORT_H_
#define _XFDTD_PORT_H_

#include <complex>
#include <memory>

#include "monitor/current_monitor.h"
#include "monitor/voltage_monitor.h"

namespace xfdtd {
class Port {
 public:
  Port(int index, std::complex<double> impedance, bool is_activated_source,
       std::shared_ptr<xfdtd::VoltageMonitor> voltage,
       std::shared_ptr<xfdtd::CurrentMonitor> current);
  Port(const Port &) = default;
  Port(Port &&) = default;

  void init(const std::shared_ptr<const FDTDBasicCoff> &fdtd_basic_coff,
            const std::shared_ptr<const GridSpace> &grid_space,
            const std::shared_ptr<const EMF> &emf);

  void calculateSParameter(const xt::xarray<double> &_frequencies);

  int getPortIndex() const;

  bool isActivatedSource() const;

  xt::xarray<std::complex<double>> getA() const;

  xt::xarray<std::complex<double>> getB() const;

 private:
  int _index;
  std::complex<double> _impedance;
  bool _is_activated_source;
  std::shared_ptr<xfdtd::VoltageMonitor> _voltage;
  std::shared_ptr<xfdtd::CurrentMonitor> _current;

  double _dt;

  xt::xarray<std::complex<double>> _a;
  xt::xarray<std::complex<double>> _b;
};

}  // namespace xfdtd

#endif  // _XFDTD_PORT_H_
