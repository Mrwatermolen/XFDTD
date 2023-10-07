#include "grid/grid.h"

namespace xfdtd {
Grid::Grid(size_t i, size_t j, size_t k) : _i{i}, _j{j}, _k{k} {}

size_t Grid::getIndexI() const { return _i; }

size_t Grid::getIndexJ() const { return _j; }

size_t Grid::getIndexK() const { return _k; }

}  // namespace xfdtd