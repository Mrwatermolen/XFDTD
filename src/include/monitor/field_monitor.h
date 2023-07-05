#ifndef _FIELD_MONITOR_H_
#define _FIELD_MONITOR_H_

#include <filesystem>
#include <memory>
#include <vector>

#include "electromagnetic_field/electromagnetic_field.h"
#include "monitor/monitor.h"
#include "shape/cube.h"
#include "shape/shape.h"
#include "util/type_define.h"

namespace xfdtd {
class TimeDomainFieldMonitor : public Monitor {
 public:
  TimeDomainFieldMonitor(std::unique_ptr<Shape> shape, PlaneType plane_type,
                         EMComponent component,
                         std::filesystem::path output_dir_path,
                         std::string _output_file_name);
  TimeDomainFieldMonitor(const TimeDomainFieldMonitor& other) = delete;
  TimeDomainFieldMonitor(TimeDomainFieldMonitor&& other) noexcept = default;
  TimeDomainFieldMonitor& operator=(const TimeDomainFieldMonitor& other) =
      delete;
  TimeDomainFieldMonitor& operator=(TimeDomainFieldMonitor&& other) noexcept =
      default;
  ~TimeDomainFieldMonitor() override = default;

  PlaneType getPlaneType() const { return _plane_type; }
  double& getEMComponent(SpatialIndex i, SpatialIndex j, SpatialIndex k);
  const YeeCellArray& getYeeCells() const { return _yee_cells; }

  void setYeeCells(const YeeCellArray& yee_cells) override;
  void setYeeCells(YeeCellArray&& yee_cells) override;

  void update(size_t current_time_step) override;

  // TODO(franzero): there are some errors when PlaneType is YZ.
  void outputData() override;

  inline void setPlaneType(PlaneType plane_type) { _plane_type = plane_type; }

  inline void setComponent(EMComponent component) { _component = component; }

 private:
  PlaneType _plane_type;
  EMComponent _component;
  YeeCellArray _yee_cells;
};
}  // namespace xfdtd

#endif  // _FIELD_MONITOR_H_
