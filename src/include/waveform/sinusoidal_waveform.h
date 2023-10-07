#ifndef _XFDTD_SINUSOIDAL_WAVEFORM_H_
#define _XFDTD_SINUSOIDAL_WAVEFORM_H_

#include "waveform/waveform.h"
namespace xfdtd {
class SinusoidalWaveform : public Waveform {
 public:
  SinusoidalWaveform(double amplitude, double frequency, double phase);
  SinusoidalWaveform(double frequency, double phase);
  SinusoidalWaveform(const SinusoidalWaveform&) = default;
  SinusoidalWaveform& operator=(const SinusoidalWaveform&) = default;
  SinusoidalWaveform(SinusoidalWaveform&&) noexcept = default;
  SinusoidalWaveform& operator=(SinusoidalWaveform&&) noexcept = default;
  ~SinusoidalWaveform() override = default;

  std::unique_ptr<Waveform> clone() const override;

  void init(xt::xarray<double> time_array) override;

  double getValueByTime(double time) const override;

 private:
  double _frequency;
  double _phase;
};
}  // namespace xfdtd

#endif  // _XFDTD_SINUSOIDAL_WAVEFORM_H_
