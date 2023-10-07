#ifndef _XFDTD_CUSTOM_WAVEFORM_H_
#define _XFDTD_CUSTOM_WAVEFORM_H_

#include <functional>

#include "waveform/waveform.h"

namespace xfdtd {
class CustomWaveform : public Waveform {
 public:
  CustomWaveform(double amplitude, std::function<double(double)> f);
  ~CustomWaveform() override = default;

  std::unique_ptr<Waveform> clone() const override;

  void init(xt::xarray<double> time_array) override;

  double getValueByTime(double time) const override;

 private:
  std::function<double(double)> _f;
};
}  // namespace xfdtd

#endif  // _XFDTD_CUSTOM_WAVEFORM_H_
