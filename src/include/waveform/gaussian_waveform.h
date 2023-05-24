#include "waveform/waveform.h"

namespace xfdtd {
class GaussianWaveform : public Waveform {
 public:
  GaussianWaveform(double amplitude, double tau, double t_0);
  GaussianWaveform(const GaussianWaveform&) = default;
  GaussianWaveform& operator=(const GaussianWaveform&) = default;
  GaussianWaveform(GaussianWaveform&&) noexcept = default;
  GaussianWaveform& operator=(GaussianWaveform&&) noexcept = default;
  ~GaussianWaveform() override = default;

  double getValueByTime(double time) const override;
  void init(const std::vector<double>& time_array) override;

 private:
  double _amplitude;
  double _tau;
  double _t_0;
};
}  // namespace xfdtd
