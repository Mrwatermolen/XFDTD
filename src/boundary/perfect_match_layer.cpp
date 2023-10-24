#include "boundary/perfect_match_layer.h"

#include <cmath>
#include <cstddef>
#include <utility>

#include "util/constant.h"
#include "util/type_define.h"

namespace xfdtd {

PML::PML(Orientation orientation, int thickness, int order, double sigma_ratio,
         double alpha_min, double alpha_max, double kappa_max)
    : _orientation{orientation},
      _thickness{thickness},
      _order{order},
      _sigma_ratio{sigma_ratio},
      _alpha_min{alpha_min},
      _alpha_max{alpha_max},
      _kappa_max{kappa_max} {}

void PML::init(double dl, double dt, SpatialIndex start_index, int na, int nb,
               EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea) {
  _dl = dl;
  _dt = dt;
  _na = na;
  _nb = nb;
  _start_index = start_index;
  _rho_e.resize({static_cast<size_t>(_thickness)});
  _rho_m.resize({static_cast<size_t>(_thickness)});
  _sigma_e.resize({static_cast<size_t>(_thickness)});
  _sigma_m.resize({static_cast<size_t>(_thickness)});
  _alpha_e.resize({static_cast<size_t>(_thickness)});
  _alpha_m.resize({static_cast<size_t>(_thickness)});
  _kappa_e.resize({static_cast<size_t>(_thickness)});
  _kappa_m.resize({static_cast<size_t>(_thickness)});
  _cpml_a_e.resize({static_cast<size_t>(_thickness)});
  _cpml_b_e.resize({static_cast<size_t>(_thickness)});
  _cpml_a_m.resize({static_cast<size_t>(_thickness)});
  _cpml_b_m.resize({static_cast<size_t>(_thickness)});

  if (isOrientatingPositive()) {
    initP(ceahb, cebha, chaeb, chbea);
  } else {
    initN(ceahb, cebha, chaeb, chbea);
  }
}

void PML::init(std::shared_ptr<EMF> emf,
               std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
               std::shared_ptr<GridSpace> grid_space,
               std::unique_ptr<Shape> shape) {
  defaultInit(emf, fdtd_basic_coff, grid_space, std::move(shape));
  auto ori{getOrientation()};
  auto dt{fdtd_basic_coff->getDt()};
  auto dx{grid_space->getGridBaseSizeX()};
  auto dy{grid_space->getGridBaseSizeY()};
  auto dz{grid_space->getGridBaseSizeZ()};
  auto nx{grid_space->getGridNumX()};
  auto ny{grid_space->getGridNumY()};
  auto nz{grid_space->getGridNumZ()};
  auto& ceyhz{fdtd_basic_coff->getCeyhz()};
  auto& cezhy{fdtd_basic_coff->getCezhy()};
  auto& chyez{fdtd_basic_coff->getChyez()};
  auto& chzey{fdtd_basic_coff->getChzey()};
  auto& cezhx{fdtd_basic_coff->getCezhx()};
  auto& cexhz{fdtd_basic_coff->getCexhz()};
  auto& chzex{fdtd_basic_coff->getChzex()};
  auto& chxez{fdtd_basic_coff->getChxez()};
  auto& cexhy{fdtd_basic_coff->getCexhy()};
  auto& ceyhx{fdtd_basic_coff->getCeyhx()};
  auto& chxey{fdtd_basic_coff->getChxey()};
  auto& chyex{fdtd_basic_coff->getChyex()};
  if (ori == Orientation::XN) {
    init(dx, dt, 0, ny, nz, ceyhz, cezhy, chyez, chzey);
    return;
  }
  if (ori == Orientation::YN) {
    init(dy, dt, 0, nz, nx, cezhx, cexhz, chzex, chxez);
    return;
  }
  if (ori == Orientation::ZN) {
    init(dz, dt, 0, nx, ny, cexhy, ceyhx, chxey, chyex);
    return;
  }
  if (ori == Orientation::XP) {
    init(dx, dt, nx - getSize(), ny, nz, ceyhz, cezhy, chyez, chzey);
    return;
  }
  if (ori == Orientation::YP) {
    init(dy, dt, ny - getSize(), nz, nx, cezhx, cexhz, chzex, chxez);
  }
  if (ori == Orientation::ZP) {
    init(dz, dt, nz - getSize(), nx, ny, cexhy, ceyhx, chxey, chyex);
  }
}

void PML::updateH() {
  auto ori{getOrientation()};
  if (ori == Orientation::XN || ori == Orientation::XP) {
    updateH(getEy(), getEz(), getHy(), getHz());
    return;
  }
  if (ori == Orientation::YN || ori == Orientation::YP) {
    updateH(getEz(), getEx(), getHz(), getHx());
    return;
  }
  if (ori == Orientation::ZN || ori == Orientation::ZP) {
    updateH(getEx(), getEy(), getHx(), getHy());
  }
}

void PML::updateE() {
  auto ori{getOrientation()};
  if (ori == Orientation::XN || ori == Orientation::XP) {
    updateE(getEy(), getEz(), getHy(), getHz());
    return;
  }
  if (ori == Orientation::YN || ori == Orientation::YP) {
    updateE(getEz(), getEx(), getHz(), getHx());
    return;
  }
  if (ori == Orientation::ZN || ori == Orientation::ZP) {
    updateE(getEx(), getEy(), getHx(), getHy());
  }
}

bool PML::isOrientatingPositive() const {
  return (_orientation == Orientation::XP || _orientation == Orientation::YP ||
          _orientation == Orientation::ZP);
}

void PML::initP(EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea) {
  double e2m{constant::MU_0 / constant::EPSILON_0};
  double sigma_max = _sigma_ratio * (_order + 1) / (150.0 * constant::PI * _dl);
  for (int i = 0; i < _thickness; ++i) {
    _rho_e(i) = (i + 1 - 0.75) / static_cast<double>(_thickness);
    _rho_m(i) = (i + 1 - 0.25) / static_cast<double>(_thickness);
    _sigma_e(i) = sigma_max * std::pow(_rho_e(i), _order);
    _sigma_m(i) = sigma_max * std::pow(_rho_m(i), _order);
    _sigma_m(i) = _sigma_m(i) * e2m;
    _kappa_e(i) = 1 + (_kappa_max - 1) * std::pow(_rho_e(i), _order);
    _kappa_m(i) = 1 + (_kappa_max - 1) * std::pow(_rho_m(i), _order);
    _alpha_e(i) = _alpha_min + (_alpha_max - _alpha_min) * (1 - _rho_e(i));
    _alpha_m(i) = _alpha_min + (_alpha_max - _alpha_min) * (1 - _rho_m(i));
    _alpha_m(i) = _alpha_m(i) * e2m;
    _cpml_b_e(i) = std::exp((-_dt / constant::EPSILON_0) *
                            (_sigma_e(i) / _kappa_e(i) + _alpha_e(i)));
    _cpml_a_e(i) = (1 / _dl) * (_cpml_b_e(i) - 1) * _sigma_e(i) /
                   (_kappa_e(i) * (_sigma_e(i) + _kappa_e(i) * _alpha_e(i)));
    _cpml_b_m(i) = std::exp((-_dt / constant::MU_0) *
                            (_sigma_m(i) / _kappa_m(i) + _alpha_m(i)));
    _cpml_a_m(i) = (1 / _dl) * (_cpml_b_m(i) - 1) * _sigma_m(i) /
                   (_kappa_m(i) * (_sigma_m(i) + _kappa_m(i) * _alpha_m(i)));
  }

  if (_orientation == Orientation::XP) {
    // a:y b:z
    _psi_ea.resize({static_cast<size_t>(_thickness), static_cast<size_t>(_na),
                    static_cast<size_t>(_nb + 1)});
    _psi_eb.resize({static_cast<size_t>(_thickness),
                    static_cast<size_t>(_na + 1), static_cast<size_t>(_nb)});
    _psi_ha.resize({static_cast<size_t>(_thickness),
                    static_cast<size_t>(_na + 1), static_cast<size_t>(_nb)});
    _psi_hb.resize({static_cast<size_t>(_thickness), static_cast<size_t>(_na),
                    static_cast<size_t>(_nb + 1)});

    _c_psi_ea.resize({static_cast<size_t>(_thickness), static_cast<size_t>(_na),
                      static_cast<size_t>(_nb + 1)});
    _c_psi_eb.resize({static_cast<size_t>(_thickness),
                      static_cast<size_t>(_na + 1), static_cast<size_t>(_nb)});
    _c_psi_ha.resize({static_cast<size_t>(_thickness),
                      static_cast<size_t>(_na + 1), static_cast<size_t>(_nb)});
    _c_psi_hb.resize({static_cast<size_t>(_thickness), static_cast<size_t>(_na),
                      static_cast<size_t>(_nb + 1)});

    // _psi_ea = std::move(allocateDoubleArray3D(_thickness, _na, _nb + 1));
    // _psi_eb = std::move(allocateDoubleArray3D(_thickness, _na + 1, _nb));
    // _psi_ha = std::move(allocateDoubleArray3D(_thickness, _na + 1, _nb));
    // _psi_hb = std::move(allocateDoubleArray3D(_thickness, _na, _nb + 1));

    // _c_psi_ea = std::move(allocateDoubleArray3D(_thickness, _na, _nb + 1));
    // _c_psi_eb = std::move(allocateDoubleArray3D(_thickness, _na + 1, _nb));
    // _c_psi_ha = std::move(allocateDoubleArray3D(_thickness, _na + 1, _nb));
    // _c_psi_hb = std::move(allocateDoubleArray3D(_thickness, _na, _nb + 1));

    for (int i = 0; i < _thickness; ++i) {
      for (int j = 0; j < _na; ++j) {
        for (int k = 0; k < _nb + 1; ++k) {
          _c_psi_ea(i, j, k) = ceahb(_start_index + i, j, k) * _dl;
          _c_psi_hb(i, j, k) = chbea(_start_index + i, j, k) * _dl;
          ceahb(_start_index + i, j, k) =
              ceahb(_start_index + i, j, k) / _kappa_e(i);
          chbea(_start_index + i, j, k) =
              chbea(_start_index + i, j, k) / _kappa_m(i);
        }
      }
      for (int j = 0; j < _na + 1; ++j) {
        for (int k = 0; k < _nb; ++k) {
          _c_psi_eb(i, j, k) = cebha(_start_index + i, j, k) * _dl;
          _c_psi_ha(i, j, k) = chaeb(_start_index + i, j, k) * _dl;
          cebha(_start_index + i, j, k) =
              cebha(_start_index + i, j, k) / _kappa_e(i);
          chaeb(_start_index + i, j, k) =
              chaeb(_start_index + i, j, k) / _kappa_m(i);
        }
      }
    }
  }

  if (_orientation == Orientation::YP) {
    _psi_ea.resize({static_cast<size_t>(_nb + 1),
                    static_cast<size_t>(_thickness), static_cast<size_t>(_na)});
    _psi_eb.resize({static_cast<size_t>(_nb), static_cast<size_t>(_thickness),
                    static_cast<size_t>(_na + 1)});
    _psi_ha.resize({static_cast<size_t>(_nb), static_cast<size_t>(_thickness),
                    static_cast<size_t>(_na + 1)});
    _psi_hb.resize({static_cast<size_t>(_nb + 1),
                    static_cast<size_t>(_thickness), static_cast<size_t>(_na)});

    _c_psi_ea.resize({static_cast<size_t>(_nb + 1),
                      static_cast<size_t>(_thickness),
                      static_cast<size_t>(_na)});
    _c_psi_eb.resize({static_cast<size_t>(_nb), static_cast<size_t>(_thickness),
                      static_cast<size_t>(_na + 1)});
    _c_psi_ha.resize({static_cast<size_t>(_nb), static_cast<size_t>(_thickness),
                      static_cast<size_t>(_na + 1)});
    _c_psi_hb.resize({static_cast<size_t>(_nb + 1),
                      static_cast<size_t>(_thickness),
                      static_cast<size_t>(_na)});

    // _psi_ea = std::move(allocateDoubleArray3D(_nb + 1, _thickness, _na));
    // _psi_eb = std::move(allocateDoubleArray3D(_nb, _thickness, _na + 1));
    // _psi_ha = std::move(allocateDoubleArray3D(_nb, _thickness, _na + 1));
    // _psi_hb = std::move(allocateDoubleArray3D(_nb + 1, _thickness, _na));

    // _c_psi_ea = std::move(allocateDoubleArray3D(_nb + 1, _thickness, _na));
    // _c_psi_eb = std::move(allocateDoubleArray3D(_nb, _thickness, _na + 1));
    // _c_psi_ha = std::move(allocateDoubleArray3D(_nb, _thickness, _na + 1));
    // _c_psi_hb = std::move(allocateDoubleArray3D(_nb + 1, _thickness, _na));

    for (int j = 0; j < _thickness; ++j) {
      for (int k = 0; k < _na; ++k) {
        for (int i = 0; i < _nb + 1; ++i) {
          _c_psi_ea(i, j, k) = ceahb(i, _start_index + j, k) * _dl;
          _c_psi_hb(i, j, k) = chbea(i, _start_index + j, k) * _dl;
          ceahb(i, _start_index + j, k) =
              ceahb(i, _start_index + j, k) / _kappa_e(j);
          chbea(i, _start_index + j, k) =
              chbea(i, _start_index + j, k) / _kappa_m(j);
        }
      }
      for (int k = 0; k < _na + 1; k++) {
        for (int i = 0; i < _nb; ++i) {
          _c_psi_eb(i, j, k) = cebha(i, _start_index + j, k) * _dl;
          _c_psi_ha(i, j, k) = chaeb(i, _start_index + j, k) * _dl;
          cebha(i, _start_index + j, k) =
              cebha(i, _start_index + j, k) / _kappa_e(j);
          chaeb(i, _start_index + j, k) =
              chaeb(i, _start_index + j, k) / _kappa_m(j);
        }
      }
    }
  }

  if (_orientation == Orientation::ZP) {
    _psi_ea = std::move(allocateDoubleArray3D(_na, _nb + 1, _thickness));
    _psi_eb = std::move(allocateDoubleArray3D(_na + 1, _nb, _thickness));
    _psi_ha = std::move(allocateDoubleArray3D(_na + 1, _nb, _thickness));
    _psi_hb = std::move(allocateDoubleArray3D(_na, _nb + 1, _thickness));

    _c_psi_ea = std::move(allocateDoubleArray3D(_na, _nb + 1, _thickness));
    _c_psi_eb = std::move(allocateDoubleArray3D(_na + 1, _nb, _thickness));
    _c_psi_ha = std::move(allocateDoubleArray3D(_na + 1, _nb, _thickness));
    _c_psi_hb = std::move(allocateDoubleArray3D(_na, _nb + 1, _thickness));

    for (int k = 0; k < _thickness; ++k) {
      for (int i = 0; i < _na; ++i) {
        for (int j = 0; j < _nb + 1; ++j) {
          _c_psi_ea(i, j, k) = ceahb(i, j, _start_index + k) * _dl;
          _c_psi_hb(i, j, k) = chbea(i, j, _start_index + k) * _dl;
          ceahb(i, j, _start_index + k) =
              ceahb(i, j, _start_index + k) / _kappa_e(k);
          chbea(i, j, _start_index + k) =
              chbea(i, j, _start_index + k) / _kappa_m(k);
        }
      }
      for (int i = 0; i < _na + 1; ++i) {
        for (int j = 0; j < _nb; ++j) {
          _c_psi_eb(i, j, k) = cebha(i, j, _start_index + k) * _dl;
          _c_psi_ha(i, j, k) = chaeb(i, j, _start_index + k) * _dl;
          cebha(i, j, _start_index + k) =
              cebha(i, j, _start_index + k) / _kappa_e(k);
          chaeb(i, j, _start_index + k) =
              chaeb(i, j, _start_index + k) / _kappa_m(k);
        }
      }
    }
  }
}

void PML::initN(EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea) {
  double e2m{constant::MU_0 / constant::EPSILON_0};
  double sigma_max = _sigma_ratio * (_order + 1) / (150.0 * constant::PI * _dl);

  for (int i = 0; i < _thickness; ++i) {
    _rho_e(i) = (_thickness - i - 0.75) / static_cast<double>(_thickness);
    _rho_m(i) = (_thickness - i - 0.25) / static_cast<double>(_thickness);
    _sigma_e(i) = sigma_max * std::pow(_rho_e(i), _order);
    _sigma_m(i) = sigma_max * std::pow(_rho_m(i), _order);
    _sigma_m(i) = _sigma_m(i) * e2m;
    _kappa_e(i) = 1 + (_kappa_max - 1) * std::pow(_rho_e(i), _order);
    _kappa_m(i) = 1 + (_kappa_max - 1) * std::pow(_rho_m(i), _order);
    _alpha_e(i) = _alpha_min + (_alpha_max - _alpha_min) * (1 - _rho_e(i));
    _alpha_m(i) = _alpha_min + (_alpha_max - _alpha_min) * (1 - _rho_m(i));
    _alpha_m(i) = _alpha_m(i) * e2m;
    _cpml_b_e(i) = std::exp((-_dt / constant::EPSILON_0) *
                            (_sigma_e(i) / _kappa_e(i) + _alpha_e(i)));
    _cpml_a_e(i) = (1 / _dl) * (_cpml_b_e(i) - 1) * _sigma_e(i) /
                   (_kappa_e(i) * (_sigma_e(i) + _kappa_e(i) * _alpha_e(i)));
    _cpml_b_m(i) = std::exp((-_dt / constant::MU_0) *
                            (_sigma_m(i) / _kappa_m(i) + _alpha_m(i)));
    _cpml_a_m(i) = (1 / _dl) * (_cpml_b_m(i) - 1) * _sigma_m(i) /
                   (_kappa_m(i) * (_sigma_m(i) + _kappa_m(i) * _alpha_m(i)));
  }

  if (_orientation == Orientation::XN) {
    // a:y b:z
    _psi_ea = std::move(allocateDoubleArray3D(_thickness, _na, _nb + 1));
    _psi_eb = std::move(allocateDoubleArray3D(_thickness, _na + 1, _nb));
    _psi_ha = std::move(allocateDoubleArray3D(_thickness, _na + 1, _nb));
    _psi_hb = std::move(allocateDoubleArray3D(_thickness, _na, _nb + 1));

    _c_psi_ea = std::move(allocateDoubleArray3D(_thickness, _na, _nb + 1));
    _c_psi_eb = std::move(allocateDoubleArray3D(_thickness, _na + 1, _nb));
    _c_psi_ha = std::move(allocateDoubleArray3D(_thickness, _na + 1, _nb));
    _c_psi_hb = std::move(allocateDoubleArray3D(_thickness, _na, _nb + 1));

    for (int i = 0; i < _thickness; ++i) {
      for (int j = 0; j < _na; ++j) {
        for (int k = 0; k < _nb + 1; ++k) {
          _c_psi_ea(i, j, k) = ceahb(_start_index + i + 1, j, k) * _dl;
          _c_psi_hb(i, j, k) = chbea(_start_index + i, j, k) * _dl;
          ceahb(_start_index + i + 1, j, k) =
              ceahb(_start_index + i + 1, j, k) / _kappa_e(i);
          chbea(_start_index + i, j, k) =
              chbea(_start_index + i, j, k) / _kappa_m(i);
        }
      }
      for (int j = 0; j < _na + 1; ++j) {
        for (int k = 0; k < _nb; ++k) {
          _c_psi_eb(i, j, k) = cebha(_start_index + i + 1, j, k) * _dl;
          _c_psi_ha(i, j, k) = chaeb(_start_index + i, j, k) * _dl;
          cebha(_start_index + i + 1, j, k) =
              cebha(_start_index + i + 1, j, k) / _kappa_e(i);
          chaeb(_start_index + i, j, k) =
              chaeb(_start_index + i, j, k) / _kappa_m(i);
        }
      }
    }
  }

  if (_orientation == Orientation::YN) {
    _psi_ea = std::move(allocateDoubleArray3D(_nb + 1, _thickness, _na));
    _psi_eb = std::move(allocateDoubleArray3D(_nb, _thickness, _na + 1));
    _psi_ha = std::move(allocateDoubleArray3D(_nb, _thickness, _na + 1));
    _psi_hb = std::move(allocateDoubleArray3D(_nb + 1, _thickness, _na));

    _c_psi_ea = std::move(allocateDoubleArray3D(_nb + 1, _thickness, _na));
    _c_psi_eb = std::move(allocateDoubleArray3D(_nb, _thickness, _na + 1));
    _c_psi_ha = std::move(allocateDoubleArray3D(_nb, _thickness, _na + 1));
    _c_psi_hb = std::move(allocateDoubleArray3D(_nb + 1, _thickness, _na));

    for (int j = 0; j < _thickness; ++j) {
      for (int k = 0; k < _na; ++k) {
        for (int i = 0; i < _nb + 1; ++i) {
          _c_psi_ea(i, j, k) = ceahb(i, _start_index + j + 1, k) * _dl;
          _c_psi_hb(i, j, k) = chbea(i, _start_index + j, k) * _dl;
          ceahb(i, _start_index + j + 1, k) =
              ceahb(i, _start_index + j + 1, k) / _kappa_e(j);
          chbea(i, _start_index + j, k) =
              chbea(i, _start_index + j, k) / _kappa_m(j);
        }
      }
      for (int k = 0; k < _na + 1; ++k) {
        for (int i = 0; i < _nb; ++i) {
          _c_psi_eb(i, j, k) = cebha(i, _start_index + j + 1, k) * _dl;
          _c_psi_ha(i, j, k) = chaeb(i, _start_index + j, k) * _dl;
          cebha(i, _start_index + j + 1, k) =
              cebha(i, _start_index + j + 1, k) / _kappa_e(j);
          chaeb(i, _start_index + j, k) =
              chaeb(i, _start_index + j, k) / _kappa_m(j);
        }
      }
    }
  }

  if (_orientation == Orientation::ZN) {
    _psi_ea = std::move(allocateDoubleArray3D(_na, _nb + 1, _thickness));
    _psi_eb = std::move(allocateDoubleArray3D(_na + 1, _nb, _thickness));
    _psi_ha = std::move(allocateDoubleArray3D(_na + 1, _nb, _thickness));
    _psi_hb = std::move(allocateDoubleArray3D(_na, _nb + 1, _thickness));

    _c_psi_ea = std::move(allocateDoubleArray3D(_na, _nb + 1, _thickness));
    _c_psi_eb = std::move(allocateDoubleArray3D(_na + 1, _nb, _thickness));
    _c_psi_ha = std::move(allocateDoubleArray3D(_na + 1, _nb, _thickness));
    _c_psi_hb = std::move(allocateDoubleArray3D(_na, _nb + 1, _thickness));

    for (int k = 0; k < _thickness; ++k) {
      for (int i = 0; i < _na; ++i) {
        for (int j = 0; j < _nb + 1; ++j) {
          _c_psi_ea(i, j, k) = ceahb(i, j, _start_index + k + 1) * _dl;
          _c_psi_hb(i, j, k) = chbea(i, j, _start_index + k) * _dl;
          ceahb(i, j, _start_index + k + 1) =
              ceahb(i, j, _start_index + k + 1) / _kappa_e(k);
          chbea(i, j, _start_index + k) =
              chbea(i, j, _start_index + k) / _kappa_m(k);
        }
      }
      for (int i = 0; i < _na + 1; ++i) {
        for (int j = 0; j < _nb; ++j) {
          _c_psi_eb(i, j, k) = cebha(i, j, _start_index + k + 1) * _dl;
          _c_psi_ha(i, j, k) = chaeb(i, j, _start_index + k) * _dl;
          cebha(i, j, _start_index + k + 1) =
              cebha(i, j, _start_index + k + 1) / _kappa_e(k);
          chaeb(i, j, _start_index + k) =
              chaeb(i, j, _start_index + k) / _kappa_m(k);
        }
      }
    }
  }
}

void PML::updateH(EFTA& ea, EFTA& eb, EFTA& ha, EFTA& hb) {
  if (_orientation == Orientation::XN) {
    for (int i = 0; i < _thickness; ++i) {
      for (int j = 0; j < _na; ++j) {
        for (int k = 0; k < _nb + 1; ++k) {
          _psi_hb(i, j, k) = _cpml_b_m(i) * _psi_hb(i, j, k) +
                             _cpml_a_m(i) * (ea(_start_index + i + 1, j, k) -
                                             ea(_start_index + i, j, k));
          hb(_start_index + i, j, k) = hb(_start_index + i, j, k) +
                                       _c_psi_hb(i, j, k) * _psi_hb(i, j, k);
        }
      }
      for (int j = 0; j < _na + 1; ++j) {
        for (int k = 0; k < _nb; ++k) {
          _psi_ha(i, j, k) = _cpml_b_m(i) * _psi_ha(i, j, k) +
                             _cpml_a_m(i) * (eb(_start_index + i + 1, j, k) -
                                             eb(_start_index + i, j, k));
          ha(_start_index + i, j, k) = ha(_start_index + i, j, k) +
                                       _c_psi_ha(i, j, k) * _psi_ha(i, j, k);
        }
      }
    }
  }

  if (_orientation == Orientation::YN) {
    for (int j = 0; j < _thickness; ++j) {
      for (int k = 0; k < _na; ++k) {
        for (int i = 0; i < _nb + 1; ++i) {
          _psi_hb(i, j, k) = _cpml_b_m(j) * _psi_hb(i, j, k) +
                             _cpml_a_m(j) * (ea(i, _start_index + j + 1, k) -
                                             ea(i, _start_index + j, k));
          hb(i, _start_index + j, k) = hb(i, _start_index + j, k) +
                                       _c_psi_hb(i, j, k) * _psi_hb(i, j, k);
        }
      }
      for (int k = 0; k < _na + 1; ++k) {
        for (int i = 0; i < _nb; ++i) {
          _psi_ha(i, j, k) = _cpml_b_m(j) * _psi_ha(i, j, k) +
                             _cpml_a_m(j) * (eb(i, _start_index + j + 1, k) -
                                             eb(i, _start_index + j, k));
          ha(i, _start_index + j, k) = ha(i, _start_index + j, k) +
                                       _c_psi_ha(i, j, k) * _psi_ha(i, j, k);
        }
      }
    }
  }

  if (_orientation == Orientation::ZN) {
    for (int k = 0; k < _thickness; ++k) {
      for (int i = 0; i < _na; ++i) {
        for (int j = 0; j < _nb + 1; ++j) {
          _psi_hb(i, j, k) = _cpml_b_m(k) * _psi_hb(i, j, k) +
                             _cpml_a_m(k) * (ea(i, j, _start_index + k + 1) -
                                             ea(i, j, _start_index + k));
          hb(i, j, _start_index + k) = hb(i, j, _start_index + k) +
                                       _c_psi_hb(i, j, k) * _psi_hb(i, j, k);
        }
      }
      for (int i = 0; i < _na + 1; ++i) {
        for (int j = 0; j < _nb; ++j) {
          _psi_ha(i, j, k) = _cpml_b_m(k) * _psi_ha(i, j, k) +
                             _cpml_a_m(k) * (eb(i, j, _start_index + k + 1) -
                                             eb(i, j, _start_index + k));
          ha(i, j, _start_index + k) = ha(i, j, _start_index + k) +
                                       _c_psi_ha(i, j, k) * _psi_ha(i, j, k);
        }
      }
    }
  }

  if (_orientation == Orientation::XP) {
    for (int i = 0; i < _thickness; ++i) {
      for (int j = 0; j < _na; ++j) {
        for (int k = 0; k < _nb + 1; ++k) {
          _psi_hb(i, j, k) = _cpml_b_m(i) * _psi_hb(i, j, k) +
                             _cpml_a_m(i) * (ea(_start_index + i + 1, j, k) -
                                             ea(_start_index + i, j, k));
          hb(_start_index + i, j, k) = hb(_start_index + i, j, k) +
                                       _c_psi_hb(i, j, k) * _psi_hb(i, j, k);
        }
      }
      for (int j = 0; j < _na + 1; ++j) {
        for (int k = 0; k < _nb; ++k) {
          _psi_ha(i, j, k) = _cpml_b_m(i) * _psi_ha(i, j, k) +
                             _cpml_a_m(i) * (eb(_start_index + i + 1, j, k) -
                                             eb(_start_index + i, j, k));
          ha(_start_index + i, j, k) = ha(_start_index + i, j, k) +
                                       _c_psi_ha(i, j, k) * _psi_ha(i, j, k);
        }
      }
    }
  }

  if (_orientation == Orientation::YP) {
    for (int j = 0; j < _thickness; ++j) {
      for (int k = 0; k < _na; ++k) {
        for (int i = 0; i < _nb + 1; ++i) {
          _psi_hb(i, j, k) = _cpml_b_m(j) * _psi_hb(i, j, k) +
                             _cpml_a_m(j) * (ea(i, _start_index + j + 1, k) -
                                             ea(i, _start_index + j, k));
          hb(i, _start_index + j, k) = hb(i, _start_index + j, k) +
                                       _c_psi_hb(i, j, k) * _psi_hb(i, j, k);
        }
      }
      for (int k = 0; k < _na + 1; ++k) {
        for (int i = 0; i < _nb; ++i) {
          _psi_ha(i, j, k) = _cpml_b_m(j) * _psi_ha(i, j, k) +
                             _cpml_a_m(j) * (eb(i, _start_index + j + 1, k) -
                                             eb(i, _start_index + j, k));
          ha(i, _start_index + j, k) = ha(i, _start_index + j, k) +
                                       _c_psi_ha(i, j, k) * _psi_ha(i, j, k);
        }
      }
    }
  }

  if (_orientation == Orientation::ZP) {
    for (int k = 0; k < _thickness; ++k) {
      for (int i = 0; i < _na; ++i) {
        for (int j = 0; j < _nb + 1; ++j) {
          _psi_hb(i, j, k) = _cpml_b_m(k) * _psi_hb(i, j, k) +
                             _cpml_a_m(k) * (ea(i, j, _start_index + k + 1) -
                                             ea(i, j, _start_index + k));
          hb(i, j, _start_index + k) = hb(i, j, _start_index + k) +
                                       _c_psi_hb(i, j, k) * _psi_hb(i, j, k);
        }
      }
      for (int i = 0; i < _na + 1; ++i) {
        for (int j = 0; j < _nb; ++j) {
          _psi_ha(i, j, k) = _cpml_b_m(k) * _psi_ha(i, j, k) +
                             _cpml_a_m(k) * (eb(i, j, _start_index + k + 1) -
                                             eb(i, j, _start_index + k));
          ha(i, j, _start_index + k) = ha(i, j, _start_index + k) +
                                       _c_psi_ha(i, j, k) * _psi_ha(i, j, k);
        }
      }
    }
  }
}

void PML::updateE(EFTA& ea, EFTA& eb, EFTA& ha, EFTA& hb) {
  if (_orientation == Orientation::XN) {
    for (int i = 0; i < _thickness; ++i) {
      for (int j = 0; j < _na; ++j) {
        for (int k = 0; k < _nb + 1; ++k) {
          _psi_ea(i, j, k) = _cpml_b_e(i) * _psi_ea(i, j, k) +
                             _cpml_a_e(i) * (hb(_start_index + i + 1, j, k) -
                                             hb(_start_index + i, j, k));
          ea(_start_index + i + 1, j, k) =
              ea(_start_index + i + 1, j, k) +
              _c_psi_ea(i, j, k) * _psi_ea(i, j, k);
        }
      }
      for (int j = 0; j < _na + 1; ++j) {
        for (int k = 0; k < _nb; ++k) {
          _psi_eb(i, j, k) = _cpml_b_e(i) * _psi_eb(i, j, k) +
                             _cpml_a_e(i) * (ha(_start_index + i + 1, j, k) -
                                             ha(_start_index + i, j, k));
          eb(_start_index + i + 1, j, k) =
              eb(_start_index + i + 1, j, k) +
              _c_psi_eb(i, j, k) * _psi_eb(i, j, k);
        }
      }
    }
  }

  if (_orientation == Orientation::YN) {
    for (int j = 0; j < _thickness; ++j) {
      for (int k = 0; k < _na; ++k) {
        for (int i = 0; i < _nb + 1; ++i) {
          _psi_ea(i, j, k) = _cpml_b_e(j) * _psi_ea(i, j, k) +
                             _cpml_a_e(j) * (hb(i, _start_index + j + 1, k) -
                                             hb(i, _start_index + j, k));
          ea(i, _start_index + j + 1, k) =
              ea(i, _start_index + j + 1, k) +
              _c_psi_ea(i, j, k) * _psi_ea(i, j, k);
        }
      }
      for (int k = 0; k < _na + 1; ++k) {
        for (int i = 0; i < _nb; ++i) {
          _psi_eb(i, j, k) = _cpml_b_e(j) * _psi_eb(i, j, k) +
                             _cpml_a_e(j) * (ha(i, _start_index + j + 1, k) -
                                             ha(i, _start_index + j, k));
          eb(i, _start_index + j + 1, k) =
              eb(i, _start_index + j + 1, k) +
              _c_psi_eb(i, j, k) * _psi_eb(i, j, k);
        }
      }
    }
  }

  if (_orientation == Orientation::ZN) {
    for (int k = 0; k < _thickness; ++k) {
      for (int i = 0; i < _na; ++i) {
        for (int j = 0; j < _nb + 1; ++j) {
          _psi_ea(i, j, k) = _cpml_b_e(k) * _psi_ea(i, j, k) +
                             _cpml_a_e(k) * (hb(i, j, _start_index + k + 1) -
                                             hb(i, j, _start_index + k));
          ea(i, j, _start_index + k + 1) =
              ea(i, j, _start_index + k + 1) +
              _c_psi_ea(i, j, k) * _psi_ea(i, j, k);
        }
      }
      for (int i = 0; i < _na + 1; ++i) {
        for (int j = 0; j < _nb; ++j) {
          _psi_eb(i, j, k) = _cpml_b_e(k) * _psi_eb(i, j, k) +
                             _cpml_a_e(k) * (ha(i, j, _start_index + k + 1) -
                                             ha(i, j, _start_index + k));
          eb(i, j, _start_index + k + 1) =
              eb(i, j, _start_index + k + 1) +
              _c_psi_eb(i, j, k) * _psi_eb(i, j, k);
        }
      }
    }
  }

  if (_orientation == Orientation::XP) {
    for (int i = 0; i < _thickness; ++i) {
      for (int j = 0; j < _na; ++j) {
        for (int k = 0; k < _nb + 1; ++k) {
          _psi_ea(i, j, k) = _cpml_b_e(i) * _psi_ea(i, j, k) +
                             _cpml_a_e(i) * (hb(_start_index + i, j, k) -
                                             hb(_start_index + i - 1, j, k));
          ea(_start_index + i, j, k) = ea(_start_index + i, j, k) +
                                       _c_psi_ea(i, j, k) * _psi_ea(i, j, k);
        }
      }
      for (int j = 0; j < _na + 1; ++j) {
        for (int k = 0; k < _nb; ++k) {
          _psi_eb(i, j, k) = _cpml_b_e(i) * _psi_eb(i, j, k) +
                             _cpml_a_e(i) * (ha(_start_index + i, j, k) -
                                             ha(_start_index + i - 1, j, k));
          eb(_start_index + i, j, k) = eb(_start_index + i, j, k) +
                                       _c_psi_eb(i, j, k) * _psi_eb(i, j, k);
        }
      }
    }
  }

  if (_orientation == Orientation::YP) {
    for (int j = 0; j < _thickness; ++j) {
      for (int k = 0; k < _na; ++k) {
        for (int i = 0; i < _nb + 1; ++i) {
          _psi_ea(i, j, k) = _cpml_b_e(j) * _psi_ea(i, j, k) +
                             _cpml_a_e(j) * (hb(i, _start_index + j, k) -
                                             hb(i, _start_index + j - 1, k));
          ea(i, _start_index + j, k) = ea(i, _start_index + j, k) +
                                       _c_psi_ea(i, j, k) * _psi_ea(i, j, k);
        }
      }
      for (int k = 0; k < _na + 1; ++k) {
        for (int i = 0; i < _nb; ++i) {
          _psi_eb(i, j, k) = _cpml_b_e(j) * _psi_eb(i, j, k) +
                             _cpml_a_e(j) * (ha(i, _start_index + j, k) -
                                             ha(i, _start_index + j - 1, k));
          eb(i, _start_index + j, k) = eb(i, _start_index + j, k) +
                                       _c_psi_eb(i, j, k) * _psi_eb(i, j, k);
        }
      }
    }
  }

  if (_orientation == Orientation::ZP) {
    for (int k = 0; k < _thickness; ++k) {
      for (int i = 0; i < _na; ++i) {
        for (int j = 0; j < _nb + 1; ++j) {
          _psi_ea(i, j, k) = _cpml_b_e(k) * _psi_ea(i, j, k) +
                             _cpml_a_e(k) * (hb(i, j, _start_index + k) -
                                             hb(i, j, _start_index + k - 1));
          ea(i, j, _start_index + k) = ea(i, j, _start_index + k) +
                                       _c_psi_ea(i, j, k) * _psi_ea(i, j, k);
        }
      }
      for (int i = 0; i < _na + 1; ++i) {
        for (int j = 0; j < _nb; ++j) {
          _psi_eb(i, j, k) = _cpml_b_e(k) * _psi_eb(i, j, k) +
                             _cpml_a_e(k) * (ha(i, j, _start_index + k) -
                                             ha(i, j, _start_index + k - 1));
          eb(i, j, _start_index + k) = eb(i, j, _start_index + k) +
                                       _c_psi_eb(i, j, k) * _psi_eb(i, j, k);
        }
      }
    }
  }
}

}  // namespace xfdtd
