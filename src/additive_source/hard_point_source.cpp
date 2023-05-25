#include "additive_source/hard_point_source.h"

#include <utility>
#include <vector>

#include "additive_source/source.h"
#include "waveform/waveform.h"

namespace xfdtd {
HardPonitSource::HardPonitSource(std::unique_ptr<Waveform> waveform,
                                 PointVector point)
    : Source(std::move(waveform)), _point{std::move(point)} {}

void HardPonitSource::init(const std::vector<double> &time_array) {
  _waveform->init(time_array);
}

PointVector HardPonitSource::getPoint() const { return _point; }
}  // namespace xfdtd
