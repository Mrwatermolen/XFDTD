#include "additive_source/source.h"

namespace xfdtd {

Source::Source(std::unique_ptr<Waveform> waveform)
    : _waveform(std::move(waveform)) {}

Source::Source(const Source &others) : _waveform(others._waveform->clone()) {}

Source::Source(Source &&others) noexcept
    : _waveform(std::move(others._waveform)) {}

Source &Source::operator=(const Source &others) {
  if (this == &others) {
    return *this;
  }

  _waveform = others._waveform->clone();
  return *this;
}

Source &Source::operator=(Source &&others) noexcept {
  if (this == &others) {
    return *this;
  }

  _waveform = std::move(others._waveform);
  return *this;
}

}  // namespace xfdtd
