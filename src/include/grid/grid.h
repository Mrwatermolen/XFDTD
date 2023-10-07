#ifndef _XFDTD_GRID_H_
#define _XFDTD_GRID_H_

#include <cstddef>

namespace xfdtd {
class Grid {
  friend class GridSpace;

 public:
  Grid(size_t i, size_t j, size_t k);
  ~Grid() = default;

  size_t getIndexI() const;

  size_t getIndexJ() const;

  size_t getIndexK() const;

  int getMaterialIndex() const { return _material_index; }

  void setMaterialIndex(int index) { _material_index = index; }

 private:
  size_t _i, _j, _k;

  int _material_index{-1};
};

}  // namespace xfdtd

#endif  // _XFDTD_GRID_H_
