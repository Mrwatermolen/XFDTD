#include "material/dispersive_material.h"

#include "material/material.h"
#include "util/constant.h"

namespace xfdtd {

LorentzMedium::LorentzMedium(std::string_view name, double eps_s,
                             double eps_inf, double omega_p, double nv,
                             double sigma_e, double sigma_m)
    : LinearDispersiveMaterial(name, 1, 1, sigma_e, sigma_m),
      _eps_s{eps_s},
      _eps_inf{eps_inf},
      _detla_eps{_eps_s - _eps_inf},
      _omega_p{omega_p},
      _nv{nv} {}

void LorentzMedium::init(double dt, double dl,
                         const std::shared_ptr<EMF> &emf) {
  LinearDispersiveMaterial::init(dt, dl, emf);
  auto eps_0{constant::EPSILON_0};
  auto sigma_e{getElectricalConductivity()};
  auto tmep_coff_0{_nv * dt + 1};

  _coff_alpha = (2 - _omega_p * _omega_p * dt * dt) / tmep_coff_0;
  _coff_xi = (_nv * dt - 1) / tmep_coff_0;
  _coff_gamma =
      (eps_0 * _detla_eps * _omega_p * _omega_p * dt * dt) / tmep_coff_0;

  auto temp_coff_1{4 * eps_0 * _eps_inf + _coff_gamma + 2 * sigma_e * dt};
  _coff_ca = (4 * eps_0 * _eps_inf - 2 * sigma_e * dt) / (temp_coff_1);
  _coff_cb = (_coff_gamma) / (temp_coff_1);
  _coff_cc = (4 * dt) / (temp_coff_1 * dl);
}

void LorentzMedium::updateEx(int i, int j, int k) {
  auto emf{getEMF()};
  auto dl{getDl()};
  auto dt{getDt()};
  auto &hy{emf->getHy()};
  auto &hz{emf->getHz()};

  auto &ex{emf->getEx(i, j, k)};
  auto &ex_prev{emf->getExPrev(i, j, k)};
  auto &jx{emf->getJx(i, j, k)};
  auto &jx_prev{emf->getJxPrev(i, j, k)};
  auto ex_old_temp{ex};
  auto jx_old_temp{jx};

  double j_sum{0};
  j_sum += 0.5 * dl * ((1 + _coff_alpha) * jx + _coff_xi * jx_prev);
  ex = _coff_ca * ex + _coff_cb * ex_prev +
       _coff_cc * (-(hy(i, j, k) - hy(i, j, k - 1)) +
                   (hz(i, j, k) - hz(i, j - 1, k)) - j_sum);
  jx = _coff_alpha * jx + _coff_xi * jx_prev +
       _coff_gamma * (ex - ex_prev) / (2 * dt);

  ex_prev = ex_old_temp;
  jx_prev = jx_old_temp;
}

void LorentzMedium::updateEy(int i, int j, int k) {
  auto emf{getEMF()};
  auto dl{getDl()};
  auto dt{getDt()};
  auto &hx{emf->getHx()};
  auto &hz{emf->getHz()};

  auto &ey{emf->getEy(i, j, k)};
  auto &ey_prev{emf->getEyPrev(i, j, k)};
  auto &jy{emf->getJy(i, j, k)};
  auto &jy_prev{emf->getJyPrev(i, j, k)};
  auto ey_old_temp{ey};
  auto jy_old_temp{jy};

  double j_sum{0};
  j_sum += 0.5 * dl * ((1 + _coff_alpha) * jy + _coff_xi * jy_prev);
  ey = _coff_ca * ey + _coff_cb * ey_prev +
       _coff_cc * (-(hz(i, j, k) - hz(i - 1, j, k)) +
                   (hx(i, j, k) - hx(i, j, k - 1)) - j_sum);

  jy = _coff_alpha * jy + _coff_xi * jy_prev +
       _coff_gamma * (ey - ey_prev) / (2 * dt);

  ey_prev = ey_old_temp;
  jy_prev = jy_old_temp;
}

void LorentzMedium::updateEz(int i, int j, int k) {
  auto emf{getEMF()};
  auto dl{getDl()};
  auto dt{getDt()};
  auto &hx{emf->getHx()};
  auto &hy{emf->getHy()};

  auto &ez{emf->getEz(i, j, k)};
  auto &ez_prev{emf->getEzPrev(i, j, k)};
  auto &jz{emf->getJz(i, j, k)};
  auto &jz_prev{emf->getJzPrev(i, j, k)};
  auto ez_old_temp{ez};
  auto jz_old_temp{jz};

  double j_sum{0};
  j_sum += 0.5 * dl * ((1 + _coff_alpha) * jz + _coff_xi * jz_prev);
  ez = _coff_ca * ez + _coff_cb * ez_prev +
       _coff_cc * (-(hx(i, j, k) - hx(i, j - 1, k)) +
                   (hy(i, j, k) - hy(i - 1, j, k)) - j_sum);

  jz = _coff_alpha * jz + _coff_xi * jz_prev +
       _coff_gamma * (ez - ez_prev) / (2 * dt);

  ez_prev = ez_old_temp;
  jz_prev = jz_old_temp;
}

DrudeMedium::DrudeMedium(std::string_view name, double eps_inf,
                         double pole_omega, double nv, double sigma_e,
                         double sigma_m)
    : LinearDispersiveMaterial(name, 1, 1, sigma_e, sigma_m),
      _eps_inf{eps_inf},
      _omega_p{pole_omega},
      _nv{nv} {}

void DrudeMedium::init(double dt, double dl, const std::shared_ptr<EMF> &emf) {
  LinearDispersiveMaterial::init(dt, dl, emf);
  auto eps_0{constant::EPSILON_0};
  auto sigma_e{getElectricalConductivity()};

  _coff_k = (2 - _nv * dt) / (2 + _nv * dt);
  _coff_beta = (_omega_p * _omega_p * eps_0 * dt) / (2 + _nv * dt);

  auto temp_coff_1{2 * eps_0 * _eps_inf + _coff_beta * dt + sigma_e * dt};
  // TODO(franzero): check this
  _coff_ca =
      (2 * eps_0 * _eps_inf - _coff_beta * dt - sigma_e * dt) / (temp_coff_1);
  _coff_cb = 2 * dt / (temp_coff_1 * dl);
}

void DrudeMedium::updateEx(int i, int j, int k) {
  auto emf{getEMF()};
  auto dl{getDl()};
  auto dt{getDt()};
  auto &hy{emf->getHy()};
  auto &hz{emf->getHz()};

  auto &ex{emf->getEx(i, j, k)};
  auto ex_old_temp{ex};
  auto &jx{emf->getJx(i, j, k)};

  double j_sum{0};
  j_sum += 0.5 * dl * ((1 + _coff_k) * jx);
  ex = _coff_ca * ex + _coff_cb * (-(hy(i, j, k) - hy(i, j, k - 1)) +
                                   (hz(i, j, k) - hz(i, j - 1, k)) - j_sum);
  jx = _coff_k * jx + _coff_beta * (ex + ex_old_temp);
}

void DrudeMedium::updateEy(int i, int j, int k) {
  auto emf{getEMF()};
  auto dl{getDl()};
  auto dt{getDt()};
  auto &hx{emf->getHx()};
  auto &hz{emf->getHz()};

  auto &ey{emf->getEy(i, j, k)};
  auto ey_old_temp{ey};
  auto &jy{emf->getJy(i, j, k)};

  double j_sum{0};
  j_sum += 0.5 * dl * ((1 + _coff_k) * jy);
  ey = _coff_ca * ey + _coff_cb * (-(hz(i, j, k) - hz(i - 1, j, k)) +
                                   (hx(i, j, k) - hx(i, j, k - 1)) - j_sum);
  jy = _coff_k * jy + _coff_beta * (ey + ey_old_temp);
}

void DrudeMedium::updateEz(int i, int j, int k) {
  auto emf{getEMF()};
  auto dl{getDl()};
  auto dt{getDt()};
  auto &hx{emf->getHx()};
  auto &hy{emf->getHy()};

  auto &ez{emf->getEz(i, j, k)};
  auto ez_old_temp{ez};
  auto &jz{emf->getJz(i, j, k)};

  double j_sum{0};
  j_sum += 0.5 * dl * ((1 + _coff_k) * jz);
  ez = _coff_ca * ez + _coff_cb * (-(hx(i, j, k) - hx(i, j - 1, k)) +
                                   (hy(i, j, k) - hy(i - 1, j, k)) - j_sum);
  jz = _coff_k * jz + _coff_beta * (ez + ez_old_temp);
}

DebyMedium::DebyMedium(std::string_view name, double eps_s, double eps_inf,
                       double tau, double sigma_e, double sigma_m)
    : LinearDispersiveMaterial(name, 1, 1, sigma_e, sigma_m),
      _eps_s{eps_s},
      _eps_inf{eps_inf},
      _detla_eps{eps_s - eps_inf},
      _tau{tau} {}

void DebyMedium::init(double dt, double dl, const std::shared_ptr<EMF> &emf) {
  LinearDispersiveMaterial::init(dt, dl, emf);
  auto eps_0{constant::EPSILON_0};
  auto sigma_e{getElectricalConductivity()};

  _coff_k = (2 * _tau - dt) / (2 * _tau + dt);
  _coff_beta = (2 * _detla_eps * eps_0 * dt) / (2 * _tau + dt);

  auto tmep_coff_1{2 * eps_0 * _eps_inf + _coff_beta + sigma_e * dt};

  _coff_ca = (2 * eps_0 * _eps_inf + _coff_beta - sigma_e * dt) / (tmep_coff_1);

  _coff_cb = (2 * dt) / ((tmep_coff_1)*dl);
}

void DebyMedium::updateEx(int i, int j, int k) {
  auto emf{getEMF()};
  auto dl{getDl()};
  auto dt{getDt()};
  auto &hy{emf->getHy()};
  auto &hz{emf->getHz()};

  auto &ex{emf->getEx(i, j, k)};
  auto ex_old_temp{ex};
  auto &jx{emf->getJx(i, j, k)};

  double j_sum{0};
  j_sum += 0.5 * dl * ((1 + _coff_k) * jx);
  ex = _coff_ca * ex + _coff_cb * (-(hy(i, j, k) - hy(i, j, k - 1)) +
                                   (hz(i, j, k) - hz(i, j - 1, k)) - j_sum);
  jx = _coff_k * jx + _coff_beta * (ex - ex_old_temp) / dt;
}

void DebyMedium::updateEy(int i, int j, int k) {
  auto emf{getEMF()};
  auto dl{getDl()};
  auto dt{getDt()};
  auto &hx{emf->getHx()};
  auto &hz{emf->getHz()};

  auto &ey{emf->getEy(i, j, k)};
  auto ey_old_temp{ey};
  auto &jy{emf->getJy(i, j, k)};

  double j_sum{0};
  j_sum += 0.5 * dl * ((1 + _coff_k) * jy);
  ey = _coff_ca * ey + _coff_cb * (-(hz(i, j, k) - hz(i - 1, j, k)) +
                                   (hx(i, j, k) - hx(i, j, k - 1)) - j_sum);
  jy = _coff_k * jy + _coff_beta * (ey - ey_old_temp) / dt;
}

void DebyMedium::updateEz(int i, int j, int k) {
  auto emf{getEMF()};
  auto dl{getDl()};
  auto dt{getDt()};
  auto &hx{emf->getHx()};
  auto &hy{emf->getHy()};

  auto &ez{emf->getEz(i, j, k)};
  auto ez_old_temp{ez};
  auto &jz{emf->getJz(i, j, k)};

  double j_sum{0};
  j_sum += 0.5 * dl * ((1 + _coff_k) * jz);
  ez = _coff_ca * ez + _coff_cb * (-(hx(i, j, k) - hx(i, j - 1, k)) +
                                   (hy(i, j, k) - hy(i - 1, j, k)) - j_sum);
  jz = _coff_k * jz + _coff_beta * (ez - ez_old_temp) / dt;
}

}  // namespace xfdtd