#ifndef _XFDTD_DISPERSIVE_MATERIAL_H_
#define _XFDTD_DISPERSIVE_MATERIAL_H_

#include <memory>
#include <xtensor/xarray.hpp>

#include "electromagnetic_field/electromagnetic_field.h"
#include "material/material.h"

namespace xfdtd {

/**
 * @brief A material's permittivity and/or permeability varies with frequency at
 * low intensities of the wave's E and H.
 *
 */
class LinearDispersiveMaterial : public Material {
 public:
  explicit LinearDispersiveMaterial(std::string_view name, float eps_r = 1,
                                    float mu_r = 1, double sigma_e = 0,
                                    double sigma_m = 0)
      : Material(name, eps_r, mu_r, sigma_e, sigma_m, true) {}

  ~LinearDispersiveMaterial() = default;

 public:
  std::shared_ptr<EMF> getEMF() { return _emf; }
  double getDl() { return _dl; }
  double getDt() { return _dt; }

  void init(double dt, double dl, const std::shared_ptr<EMF>& emf) override {
    _emf = emf;
    _dl = dl;
    _dt = dt;
  };

 private:
  std::shared_ptr<EMF> _emf;
  double _dl;
  double _dt;
};

class LorentzMedium : public LinearDispersiveMaterial {
 public:
  /**
   * @brief Construct a new Lorentz Medium object
   *
   * @param name the name of the material
   * @param eps_s the static or zero-frequncy relative permittivity
   * @param eps_inf the relative permittivity at infinite frequency
   * @param omega_p the frequency of the pole pair
   * @param nv the pole relaxation time
   * @param sigma_e the electrical conductivity
   * @param sigma_m the magnetic conductivity
   */
  LorentzMedium(std::string_view name, double eps_s, double eps_inf,
                double omega_p, double nv, double sigma_e = 0,
                double sigma_m = 0);

  ~LorentzMedium() = default;

  void init(double dt, double dl, const std::shared_ptr<EMF>& emf) override;

  void updateEx(int i, int j, int k) override;
  void updateEy(int i, int j, int k) override;
  void updateEz(int i, int j, int k) override;

  inline double getCoffAlpha() { return _coff_alpha; }
  inline double getCoffXi() { return _coff_xi; }
  inline double getCoffGamma() { return _coff_gamma; }
  inline double getCoffCa() { return _coff_ca; }
  inline double getCoffCb() { return _coff_cb; }
  inline double getCoffCc() { return _coff_cc; }

 private:
  double _eps_s;    // the static or zero-frequncy relative permittivity
  double _eps_inf;  // the relative permittivity at infinite frequency
  double _omega_p;  // the frequency of the pole pair
  double _nv;       // the pole relaxation time
  double _detla_eps;
  double _coff_alpha, _coff_xi, _coff_gamma, _coff_ca, _coff_cb, _coff_cc;
};

class DrudeMedium : public LinearDispersiveMaterial {
 public:
  /**
   * @brief Construct a new Drude Medium object
   *
   * @param pole_omega the frequency of the pole pair
   * @param nv the inverse of the pole relaxation time
   * @param sigma_e the electrical conductivity
   * @param sigma_m the magnetic conductivity
   */
  DrudeMedium(std::string_view name, double eps_inf, double pole_omega,
              double nv, double sigma_e = 0, double sigma_m = 0);

  ~DrudeMedium() = default;

  void init(double dt, double dl, const std::shared_ptr<EMF>& emf) override;

  void updateEx(int i, int j, int k) override;
  void updateEy(int i, int j, int k) override;
  void updateEz(int i, int j, int k) override;

  double getCoffK() { return _coff_k; }
  double getCoffBeta() { return _coff_beta; }
  double getCoffCa() { return _coff_ca; }
  double getCoffCb() { return _coff_cb; }

 private:
  double _eps_inf;  // the relative permittivity at infinite frequency
  double _omega_p;  // the frequency of the pole pair
  double _nv;       // the inverse of the pole relaxation time
  double _coff_k, _coff_beta, _coff_ca, _coff_cb;
};

class DebyMedium : public LinearDispersiveMaterial {
 public:
  DebyMedium(std::string_view name, double eps_s, double eps_inf, double tau,
             double sigma_e = 0, double sigma_m = 0);

  ~DebyMedium() = default;

  void init(double dt, double dl, const std::shared_ptr<EMF>& emf) override;

  void updateEx(int i, int j, int k) override;
  void updateEy(int i, int j, int k) override;
  void updateEz(int i, int j, int k) override;

  double getCoffK() { return _coff_k; }
  double getCoffBeta() { return _coff_beta; }
  double getCoffCa() { return _coff_ca; }
  double getCoffCb() { return _coff_cb; }

 private:
  double _eps_s;    // the static or zero-frequncy relative permittivity
  double _eps_inf;  // the relative permittivity at infinite frequency
  double _detla_eps;
  double _tau;  // the pole relaxation time
  double _coff_k, _coff_beta, _coff_ca, _coff_cb;
};
}  // namespace xfdtd

#endif  // _XFDTD_DISPERSIVE_MATERIAL_H_