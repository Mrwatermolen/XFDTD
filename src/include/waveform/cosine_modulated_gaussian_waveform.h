#ifndef __COSINE_MODULATED_GAUSSIAN_WAVEFORM_H__
#define __COSINE_MODULATED_GAUSSIAN_WAVEFORM_H__

#include <vector>

#include "waveform/waveform.h"
namespace xfdtd {
class CosineModulatedGaussianWaveform : public Waveform {
 public:
  CosineModulatedGaussianWaveform(double amplitude, double tau, double t_0,
                                  double modulation_frequency);
  ~CosineModulatedGaussianWaveform() override = default;
  double getValueByTime(double time) const override;
  void init(const std::vector<double>& time_array) override;

 private:
  double _amplitude;
  double _tau;
  double _t_0;
  double _modulation_frequency;
};
}  // namespace xfdtd

#endif  // __COSINE_MODULATED_GAUSSIAN_WAVEFORM_H__
