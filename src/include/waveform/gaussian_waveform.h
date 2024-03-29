#ifndef _XFDTD_GAUSSIAN_WAVEFORM_H_
#define _XFDTD_GAUSSIAN_WAVEFORM_H_

#include "waveform/waveform.h"

namespace xfdtd {
class GaussianWaveform : public Waveform {
 public:
  GaussianWaveform(double amplitude, double tau, double t_0);
  GaussianWaveform(double tau, double t_0);
  GaussianWaveform(const GaussianWaveform&) = default;
  GaussianWaveform& operator=(const GaussianWaveform&) = default;
  GaussianWaveform(GaussianWaveform&&) noexcept = default;
  GaussianWaveform& operator=(GaussianWaveform&&) noexcept = default;
  ~GaussianWaveform() override = default;

  std::unique_ptr<Waveform> clone() const override;

  void init(xt::xarray<double> time_array) override;

  double getValueByTime(double time) const override;

 private:
  double _tau;
  double _t_0;
};
}  // namespace xfdtd

#endif  // _XFDTD_GAUSSIAN_WAVEFORM_H_
