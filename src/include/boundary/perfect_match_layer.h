#ifndef _PERFECT_MATCH_LAYER_H_
#define _PERFECT_MATCH_LAYER_H_

#include <Eigen/Core>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>

#include "boundary/boundary.h"
#include "shape/shape.h"
#include "simulation/simulation.h"

namespace xfdtd {
class PML : public Boundary {
 public:
  PML(Orientation orientation, int thickness, int order = 4,
      double sigma_ratio = 1, double alpha_min = 0, double alpha_max = 0.05,
      double kappa_max = 10);
  PML(const PML& pml) = default;
  PML& operator=(const PML& pml) = default;
  PML(PML&& pml) noexcept = default;
  PML& operator=(PML&& pml) noexcept = default;
  ~PML() override = default;

  void initSize(double dl);
  void init(double dl, double dt, SpatialIndex start_index, int na, int nb,
            EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea);

  inline int getSize() const override { return _thickness; }
  inline Orientation getOrientation() const override { return _orientation; }

  const Eigen::ArrayXd& getKappaE() const;
  const Eigen::ArrayXd& getKappaM() const;
  bool isOrientatingPositive() const;

  void updateE(EFTA& ea, EFTA& eb, EFTA& ha, EFTA& hb);
  void updateH(EFTA& ea, EFTA& eb, EFTA& ha, EFTA& hb);

 private:
  Orientation _orientation;
  int _thickness;
  int _order;
  double _sigma_ratio;
  double _alpha_min;
  double _alpha_max;
  double _kappa_max;

  // init parameters
  double _dl;
  double _dt;
  SpatialIndex _na;
  SpatialIndex _nb;
  SpatialIndex _start_index;

  Eigen::ArrayXd _rho_e;
  Eigen::ArrayXd _rho_m;
  Eigen::ArrayXd _sigma_e;
  Eigen::ArrayXd _sigma_m;
  Eigen::ArrayXd _kappa_e;
  Eigen::ArrayXd _kappa_m;
  Eigen::ArrayXd _alpha_e;
  Eigen::ArrayXd _alpha_m;

  // update parameters
  Eigen::ArrayXd _cpml_a_e;
  Eigen::ArrayXd _cpml_b_e;
  Eigen::ArrayXd _cpml_a_m;
  Eigen::ArrayXd _cpml_b_m;

  Eigen::Tensor<double, 3> _psi_ea;
  Eigen::Tensor<double, 3> _psi_eb;
  Eigen::Tensor<double, 3> _psi_ha;
  Eigen::Tensor<double, 3> _psi_hb;

  Eigen::Tensor<double, 3> _c_psi_ea;
  Eigen::Tensor<double, 3> _c_psi_eb;
  Eigen::Tensor<double, 3> _c_psi_ha;
  Eigen::Tensor<double, 3> _c_psi_hb;

  void initP(EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea);
  void initN(EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea);
};
}  // namespace xfdtd
#endif  // _PERFECT_MATCH_LAYER_H_
