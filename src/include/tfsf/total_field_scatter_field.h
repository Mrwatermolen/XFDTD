#include "additive_source/source.h"
#include "shape/shape.h"

namespace xfdtd {
class TFSF : public Source {
 public:
  TFSF(std::unique_ptr<Waveform> waveform, std::unique_ptr<Shape> shape,
       std::unique_ptr<Shape> monitor);
  TFSF(const TFSF &others);
  TFSF(TFSF &&others) noexcept;
  TFSF &operator=(const TFSF &others);
  TFSF &operator=(TFSF &&others) noexcept;
  ~TFSF() override = default;

  void init(const std::vector<double> &time_array) override;

 protected:
  std::unique_ptr<Shape> _shape;
};
}  // namespace xfdtd