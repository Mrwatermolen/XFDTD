#include "waveform/waveform.h"

#include <filesystem>
#include <fstream>
#include <utility>
#include <xtensor/xcsv.hpp>
#include <xtensor/xnpy.hpp>
#include <xtensor/xview.hpp>

namespace xfdtd {

Waveform::Waveform(double amplitude) : _amplitude(amplitude) {}

double Waveform::getAmplitude() const { return _amplitude; }

void Waveform::setAmplitude(double amplitude) { _amplitude = amplitude; }

void Waveform::setValue(xt::xarray<double> value) { _value = std::move(value); }

void Waveform::dumpToCsv(const std::filesystem::path &save_path) const {
  std::ofstream file{save_path};
  auto temp_value = xt::view(_value, xt::newaxis(), xt::all());
  xt::dump_csv(file, temp_value);
  file.close();
}

void Waveform::dumpToNpy(const std::filesystem::path &save_path) const {
  xt::dump_npy(save_path, _value);
}

// void Waveform::dumpToNpy(std::string_view file_name) const {

// }
}  // namespace xfdtd
