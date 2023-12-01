#ifndef _XFDTD_FIELD_MONITOR_H_
#define _XFDTD_FIELD_MONITOR_H_

#include <memory>

#include "electromagnetic_field/electromagnetic_field.h"
#include "grid/grid_box.h"
#include "monitor/monitor.h"
#include "shape/cube.h"
#include "util/type_define.h"

namespace xfdtd {

class TimeDomainFieldMonitor : public Monitor {
 public:
  TimeDomainFieldMonitor(std::unique_ptr<Cube> shape, PlaneType plane_type,
                         EMComponent component, std::string output_dir_path,
                         std::string output_file_name);

  ~TimeDomainFieldMonitor() override = default;

  std::unique_ptr<Monitor> clone() const override;

  void init(const std::shared_ptr<const FDTDBasicCoff> &fdtd_basic_coff,
            const std::shared_ptr<const GridSpace> &grid_space,
            const std::shared_ptr<const EMF> &emf) override;

  void update() override;

  void outputData() override;

 private:
  PlaneType _plane_type;
  EMComponent _component;
  GridBox _grid_box;

  xt::xarray<double> _coord;
  xt::xarray<double> _x, _y, _z;
};

}  // namespace xfdtd

#endif  // _XFDTD_FIELD_MONITOR_H_
