#ifndef __SINUSOIDAL_WAVEFORM_H__
#define __SINUSOIDAL_WAVEFORM_H__

#include "waveform/waveform.h"
namespace xfdtd {
class SinusoidalWaveform : public Waveform {
 public:
  SinusoidalWaveform(double amplitude, double frequency, double phase);
  SinusoidalWaveform(const SinusoidalWaveform&) = default;
  SinusoidalWaveform& operator=(const SinusoidalWaveform&) = default;
  SinusoidalWaveform(SinusoidalWaveform&&) noexcept = default;
  SinusoidalWaveform& operator=(SinusoidalWaveform&&) noexcept = default;
  ~SinusoidalWaveform() override = default;

  void init(const std::vector<double>& time_array) override;
  std::unique_ptr<Waveform> clone() const override;

 private:
  double _amplitude;
  double _frequency;
  double _phase;
};
}  // namespace xfdtd

#endif
