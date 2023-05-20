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

  std::shared_ptr<double[]> getAllValues() const;
  std::shared_ptr<double[]> getValues(size_t start, size_t end) const;
  double getValue(size_t time_step) const;
  virtual double getValueByTime(double time) const = 0;

  virtual void init(const std::vector<double>& time_array) = 0;
  virtual std::unique_ptr<Waveform> clone() const = 0;

 protected:
  size_t _size;
  std::shared_ptr<double[]> _values;  // waveform values
 private:
};
}  // namespace xfdtd

#endif  // _WAVEFORM_H_
