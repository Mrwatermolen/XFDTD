#ifndef _XFDTD_GRID_BOX_H_
#define _XFDTD_GRID_BOX_H_

#include <cstddef>
#include <string>

namespace xfdtd {
class GridBox {
 public:
  GridBox() = default;

  GridBox(size_t origin_i, size_t origin_j, size_t origin_k, size_t nx,
          size_t ny, size_t nz);

  ~GridBox() = default;

  enum class Type { UNDEFINED_TYPE, POINT, LINE, FACE, CUBE };

  static std::string toString(Type dimension);

  Type getType() const;

  size_t getGridNumX() const;

  size_t getGridNumY() const;

  size_t getGridNumZ() const;

  size_t getGridNum() const;

  size_t getGridOriginIndexX() const;

  size_t getGridOriginIndexY() const;

  size_t getGridOriginIndexZ() const;

  size_t getGridStartIndexX() const;

  size_t getGridStartIndexY() const;

  size_t getGridStartIndexZ() const;

  size_t getGridEndIndexX() const;

  size_t getGridEndIndexY() const;

  size_t getGridEndIndexZ() const;

  size_t getGridCenterIndexX() const;

  size_t getGridCenterIndexY() const;

  size_t getGridCenterIndexZ() const;

 private:
  Type _dimension{Type::UNDEFINED_TYPE};
  size_t _origin_i, _origin_j, _origin_k;
  size_t _nx, _ny, _nz;
};
}  // namespace xfdtd

#endif  // _XFDTD_GRID_BOX_H_
