#ifndef __HARD_POINT_SOURCE_H__
#define __HARD_POINT_SOURCE_H__

#include <Eigen/Core>
#include <memory>

#include "additive_source/source.h"
#include "waveform/waveform.h"

namespace xfdtd {

// only for 1d test
class HardPonitSource : public Source {
 public:
  HardPonitSource(std::unique_ptr<Waveform> waveform, Eigen::Vector3d point);
  ~HardPonitSource() override = default;

  void init(const std::vector<double>& time_array) override;

  Eigen::Vector3d getPoint() const;

  double getJe() const;

 private:
  Eigen::Vector3d _point;
};
}  // namespace xfdtd

#endif
