#ifndef _XFDTD_COSINE_MODULATED_GAUSSIAN_WAVEFORM_H_
#define _XFDTD_COSINE_MODULATED_GAUSSIAN_WAVEFORM_H_

#include "waveform/waveform.h"
namespace xfdtd {
class CosineModulatedGaussianWaveform : public Waveform {
 public:
  CosineModulatedGaussianWaveform(double amplitude, double tau, double t_0,
                                  double modulation_frequency);
  ~CosineModulatedGaussianWaveform() override = default;
  double getValueByTime(double time) const override;

 private:
  double _amplitude;
  double _tau;
  double _t_0;
  double _modulation_frequency;
};
}  // namespace xfdtd

#endif  // _XFDTD_COSINE_MODULATED_GAUSSIAN_WAVEFORM_H_
