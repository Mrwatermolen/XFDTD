#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include <memory>
#include <string>
#include <string_view>
#include <tuple>

#include "electromagnetic_field/electromagnetic_field.h"

namespace xfdtd {
class Material {
 public:
  Material(std::string_view name, float eps_r, float mu_r, double sigma_e,
           double sigma_m, bool is_dispersion = false);
  Material(const Material& material) = default;
  Material(Material&& material) = default;
  Material& operator=(const Material& material) = default;
  Material& operator=(Material&& material) noexcept = default;
  ~Material() = default;

  explicit operator std::string() const;

  virtual void init(double dt, double dl, const std::shared_ptr<EMF>& emf){};

  virtual void updateEx(int i, int j, int k){};

  virtual void updateEy(int i, int j, int k){};

  virtual void updateEz(int i, int j, int k){};

  virtual void updateHx(size_t i, size_t j, size_t k){};

  virtual void updateHy(size_t i, size_t j, size_t k){};

  virtual void updateHz(size_t i, size_t j, size_t k){};

  inline double getPermittivityE() { return _eps; }
  inline double getPermeabilityM() { return _mu; }
  inline double getElectricalConductivity() { return _sigma_e; }
  inline double getMagneticConductivity() { return _sigma_m; }

  /**
   * @brief Get the Electromagnetic Properties object
   *
   * @return std::tuple<double, double, double, double> eps mu sigma_e
   * sigma_m
   */
  std::tuple<double, double, double, double> getElectromagneticProperties() {
    return {getPermittivityE(), getPermeabilityM(), getElectricalConductivity(),
            getMagneticConductivity()};
  }

  inline bool isDispersion() { return _is_dispersion; }

  static std::unique_ptr<Material> createAir(
      const std::string& name = "XFDTD_AIR");

  static std::unique_ptr<Material> createPEC();

  static std::unique_ptr<Material> createPMC();

 private:
  std::string _name;

  float _eps_r;
  float _mu_r;

  double _eps;      // permittivity
  double _mu;       // permeability
  double _sigma_e;  // electrical conductivity
  double _sigma_m;  // magnetic conductivity
  bool _is_dispersion;
};
}  // namespace xfdtd

#endif  // _MATERIAL_H_
