#ifndef __HARD_POINT_SOURCE_H__
#define __HARD_POINT_SOURCE_H__

#include <memory>

#include "additive_source/source.h"
#include "util/type_define.h"
#include "waveform/waveform.h"

namespace xfdtd {

// only for 1d test
class HardPonitSource : public Source {
 public:
  HardPonitSource(std::unique_ptr<Waveform> waveform, PointVector point);
  HardPonitSource(const HardPonitSource& other) = delete;
  HardPonitSource(HardPonitSource&& other) noexcept = default;
  ~HardPonitSource() override = default;

  void init(const std::vector<double>& time_array) override;

  PointVector getPoint() const;

  double getJe() const;

 private:
  PointVector _point;
};
}  // namespace xfdtd

#endif
