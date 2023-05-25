#ifndef _PERFECT_MATCH_LAYER_H_
#define _PERFECT_MATCH_LAYER_H_

#include <Eigen/Core>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

#include "boundary/boundary.h"
#include "shape/shape.h"
#include "util/type_define.h"

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

  void init(Simulation* simulation) override;

  inline SpatialIndex getSize() const override { return _thickness; }
  inline Orientation getOrientation() const override { return _orientation; }
  inline EFTA& getEx() { return getEMFInstance()->getEx(); }
  inline EFTA& getEy() { return getEMFInstance()->getEy(); }
  inline EFTA& getEz() { return getEMFInstance()->getEz(); }
  inline EFTA& getHx() { return getEMFInstance()->getHx(); }
  inline EFTA& getHy() { return getEMFInstance()->getHy(); }
  inline EFTA& getHz() { return getEMFInstance()->getHz(); }

  const DoubleArrary1D& getKappaE() const;
  const DoubleArrary1D& getKappaM() const;
  bool isOrientatingPositive() const;

  void updateH() override;
  void updateE() override;

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

  DoubleArrary1D _rho_e;
  DoubleArrary1D _rho_m;
  DoubleArrary1D _sigma_e;
  DoubleArrary1D _sigma_m;
  DoubleArrary1D _kappa_e;
  DoubleArrary1D _kappa_m;
  DoubleArrary1D _alpha_e;
  DoubleArrary1D _alpha_m;

  // update parameters
  DoubleArrary1D _cpml_a_e;
  DoubleArrary1D _cpml_b_e;
  DoubleArrary1D _cpml_a_m;
  DoubleArrary1D _cpml_b_m;

  DoubleArrary3D _psi_ea;
  DoubleArrary3D _psi_eb;
  DoubleArrary3D _psi_ha;
  DoubleArrary3D _psi_hb;

  DoubleArrary3D _c_psi_ea;
  DoubleArrary3D _c_psi_eb;
  DoubleArrary3D _c_psi_ha;
  DoubleArrary3D _c_psi_hb;

  void init(double dl, double dt, SpatialIndex start_index, int na, int nb,
            EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea);

  void initP(EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea);
  void initN(EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea);

  void updateE(EFTA& ea, EFTA& eb, EFTA& ha, EFTA& hb);
  void updateH(EFTA& ea, EFTA& eb, EFTA& ha, EFTA& hb);
};
}  // namespace xfdtd
#endif  // _PERFECT_MATCH_LAYER_H_
