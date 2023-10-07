#ifndef _XFDTD_WAVEFORM_H_
#define _XFDTD_WAVEFORM_H_

#include <filesystem>
#include <memory>
#include <utility>
#include <xtensor/xarray.hpp>
namespace xfdtd {
class Waveform {
 public:
  explicit Waveform(double amplitude);
  Waveform() = default;
  Waveform(const Waveform&) = default;
  Waveform& operator=(const Waveform&) = default;
  Waveform(Waveform&&) noexcept = default;
  Waveform& operator=(Waveform&&) noexcept = default;
  virtual ~Waveform() = default;

  template <typename... Args>
  double operator()(Args&&... args);

  virtual std::unique_ptr<Waveform> clone() const = 0;

  virtual void init(xt::xarray<double> time_array) = 0;

  virtual double getValueByTime(double time) const = 0;

  double getAmplitude() const;

  void setAmplitude(double amplitude);

  void dumpToCsv(const std::filesystem::path& save_path) const;

  void dumpToNpy(const std::filesystem::path& save_path) const;

 protected:
  void setValue(xt::xarray<double> value);

 private:
  double _amplitude{1.0};
  xt::xarray<double> _value;
};

template <typename... Args>
inline double Waveform::operator()(Args&&... args) {
  return _value(std::forward<Args>(args)...);
}

}  // namespace xfdtd

#endif  // _XFDTD_WAVEFORM_H_
