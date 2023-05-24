#ifndef _SOURCE_H_
#define _SOURCE_H_

#include <memory>

#include "waveform/waveform.h"
namespace xfdtd {
class Source {
 public:
  explicit Source(std::unique_ptr<Waveform> waveform);
  Source(const Source &others) = delete;
  Source(Source &&others) noexcept = default;
  Source &operator=(const Source &others) = delete;
  Source &operator=(Source &&others) noexcept = default;
  virtual ~Source() = default;

  virtual void init(const std::vector<double> &time_array) = 0;
  inline double getValue(size_t time_step) {
    return _waveform->getValue(time_step);
  }

 protected:
  std::unique_ptr<Waveform> _waveform;
};
}  // namespace xfdtd

#endif  // _SOURCE_H_
