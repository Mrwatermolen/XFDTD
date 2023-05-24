#include "additive_source/source.h"

namespace xfdtd {

Source::Source(std::unique_ptr<Waveform> waveform)
    : _waveform(std::move(waveform)) {}

}  // namespace xfdtd
