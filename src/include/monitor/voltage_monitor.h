#ifndef _XFDTD_VOLTAGE_MONITOR_H_
#define _XFDTD_VOLTAGE_MONITOR_H_

#include <cstddef>
#include <memory>

#include "monitor/monitor.h"
#include "shape/cube.h"
#include "util/type_define.h"

namespace xfdtd {

class VoltageMonitor : public Monitor {
 public:
  VoltageMonitor(std::unique_ptr<Cube> shape, Orientation orientation,
                 std::filesystem::path output_dir_path,
                 std::string output_file_name);
  ~VoltageMonitor() override = default;

  std::unique_ptr<Monitor> clone() const override;

  void init(const std::shared_ptr<const FDTDBasicCoff> &fdtd_basic_coff,
            const std::shared_ptr<const GridSpace> &grid_space,
            const std::shared_ptr<const EMF> &emf) override;

  void update() override;

  void outputData() override;

  double getValue(size_t i) const;

  auto getValueArray() const;

 private:
  Orientation _orientation;
  size_t _is, _ie, _js, _je, _ks, _ke;
  // double _dc;
  xt::xarray<double> _dc;
  xt::xarray<double> _coff;

  xt::xarray<double> _value;
};

inline auto VoltageMonitor::getValueArray() const { return _value; }

}  // namespace xfdtd

#endif  // _XFDTD_VOLTAGE_MONITOR_H_
