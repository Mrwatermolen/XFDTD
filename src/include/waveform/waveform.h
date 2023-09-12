#ifndef _XFDTD_WAVEFORM_H_
#define _XFDTD_WAVEFORM_H_

namespace xfdtd {
class Waveform {
 public:
  Waveform() = default;
  Waveform(const Waveform&) = default;
  Waveform& operator=(const Waveform&) = default;
  Waveform(Waveform&&) noexcept = default;
  Waveform& operator=(Waveform&&) noexcept = default;
  virtual ~Waveform() = default;

  virtual double getValueByTime(double time) const = 0;
};
}  // namespace xfdtd

#endif  // _XFDTD_WAVEFORM_H_
