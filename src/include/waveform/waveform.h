#ifndef _WAVEFORM_H_
#define _WAVEFORM_H_

#include <cstddef>
#include <memory>
#include <vector>
namespace xfdtd {
class Waveform {
 public:
  Waveform() = default;
  Waveform(const Waveform&) = default;
  Waveform& operator=(const Waveform&) = default;
  Waveform(Waveform&&) noexcept = default;
  Waveform& operator=(Waveform&&) noexcept = default;
  virtual ~Waveform() = default;

  const std::vector<double>& getAllValues() const;
  std::vector<double>& getAllValues();
  double getValue(size_t time_step) const;
  virtual double getValueByTime(double time) const = 0;

  virtual void init(const std::vector<double>& time_array) = 0;

 private:
  std::vector<double> _values;  // waveform values
};
}  // namespace xfdtd

#endif  // _WAVEFORM_H_
