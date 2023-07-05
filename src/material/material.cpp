#include "material/material.h"

#include "util/constant.h"

namespace xfdtd {

Material::Material(std::string_view name, float eps_r, float mu_r,
                   double sigma_e, double sigma_m, bool is_dispersion)
    : _name{name},
      _eps_r{eps_r},
      _mu_r{mu_r},
      _eps{_eps_r * constant::EPSILON_0},
      _mu{_mu_r * constant::MU_0},
      _sigma_e{sigma_e},
      _sigma_m{sigma_m},
      _is_dispersion{is_dispersion} {}

Material::operator std::string() const {
  return std::string("Material: ") + _name + "\n" +
         "eps_r: " + std::to_string(_eps_r) + "\n" +
         "mu_r: " + std::to_string(_mu_r) + "\n" +
         "eps: " + std::to_string(_eps) + "\n" + "mu: " + std::to_string(_mu) +
         "\n" + "sigma_e: " + std::to_string(_sigma_e) + "\n" +
         "sigma_m: " + std::to_string(_sigma_m) + "\n" +
         "is_dispersion: " + std::to_string(static_cast<int>(_is_dispersion));
}
}  // namespace xfdtd
