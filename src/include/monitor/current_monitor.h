#ifndef _XFDTD_CURRENT_MONITOR_H_
#define _XFDTD_CURRENT_MONITOR_H_

#include "monitor/monitor.h"

namespace xfdtd {

class CurrentMonitor : public Monitor {
 public:
  CurrentMonitor(std::unique_ptr<Cube> shape, Orientation orientation,
                 std::string output_dir_path, std::string output_file_name);
  ~CurrentMonitor() override = default;

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
  // double _da, _db;
  xt::xarray<double> _da, _db, _dc;
  double _positive;

  xt::xarray<double> _value;
};

inline auto CurrentMonitor::getValueArray() const { return _value; }

}  // namespace xfdtd

#endif  // _XFDTD_CURRENT_MONITOR_H_
