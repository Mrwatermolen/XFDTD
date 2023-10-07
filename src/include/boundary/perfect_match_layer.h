#ifndef _PERFECT_MATCH_LAYER_H_
#define _PERFECT_MATCH_LAYER_H_

#include <memory>

#include "boundary/boundary.h"
#include "shape/shape.h"
#include "util/type_define.h"

namespace xfdtd {

/**
 * @brief Perfect Match Layer.
 * This library uses the CPML(Coordinate stretched Perfectly Matched Layer)
 * method.
 * TODO(franzero): refactor the PML class.
 *
 */
class PML : public Boundary {
 public:
  /**
   * @brief Construct a new PML object
   *
   * @param orientation Normal vector of the boundary
   * @param thickness The thickness of the boundary, also called the number of
   * cells.
   * @param order The order of the CPML, default is 4.
   * @param sigma_ratio
   * @param alpha_min
   * @param alpha_max
   * @param kappa_max
   */
  PML(Orientation orientation, int thickness, int order = 4,
      double sigma_ratio = 1, double alpha_min = 0, double alpha_max = 0.05,
      double kappa_max = 10);
  PML(const PML& pml) = delete;
  PML& operator=(const PML& pml) = delete;
  PML(PML&& pml) noexcept = default;
  PML& operator=(PML&& pml) noexcept = default;
  ~PML() override = default;

  void init(std::shared_ptr<EMF> emf,
            std::shared_ptr<FDTDBasicCoff> fdtd_basic_coff,
            std::shared_ptr<GridSpace> grid_space,
            std::unique_ptr<Shape> shape) override;

  inline SpatialIndex getSize() const override { return _thickness; }
  inline Orientation getOrientation() const override { return _orientation; }

  /**
   * @brief Determine the positive and negative sign of the normal vector
   *
   * @return true The normal vector is positive(+).
   * @return false The normal vector is negative(-).
   */
  bool isOrientatingPositive() const;

  void updateH() override;
  void updateE() override;

 private:
  Orientation _orientation;  // normal vector
  int _thickness;            // number of cells
  int _order;                // order of the CPML
  double _sigma_ratio;
  double _alpha_min;
  double _alpha_max;
  double _kappa_max;

  double _dl;
  double _dt;

  // The cross product of the unit vector a and the unit vector b is the normal
  // vector of the boundary. This is the meaning of a and b.
  SpatialIndex _na;
  SpatialIndex _nb;
  SpatialIndex _start_index;

  xt::xarray<double> _rho_e;
  xt::xarray<double> _rho_m;
  xt::xarray<double> _sigma_e;
  xt::xarray<double> _sigma_m;
  xt::xarray<double> _kappa_e;
  xt::xarray<double> _kappa_m;
  xt::xarray<double> _alpha_e;
  xt::xarray<double> _alpha_m;

  // update parameters
  xt::xarray<double> _cpml_a_e;
  xt::xarray<double> _cpml_b_e;
  xt::xarray<double> _cpml_a_m;
  xt::xarray<double> _cpml_b_m;

  xt::xtensor<double, 3> _psi_ea;
  xt::xtensor<double, 3> _psi_eb;
  xt::xtensor<double, 3> _psi_ha;
  xt::xtensor<double, 3> _psi_hb;

  xt::xtensor<double, 3> _c_psi_ea;
  xt::xtensor<double, 3> _c_psi_eb;
  xt::xtensor<double, 3> _c_psi_ha;
  xt::xtensor<double, 3> _c_psi_hb;

  void init(double dl, double dt, SpatialIndex start_index, int na, int nb,
            EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea);

  // TODO(franzero): These two functions can be combined.
  void initP(EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea);
  void initN(EFTA& ceahb, EFTA& cebha, EFTA& chaeb, EFTA& chbea);

  void updateE(EFTA& ea, EFTA& eb, EFTA& ha, EFTA& hb);
  void updateH(EFTA& ea, EFTA& eb, EFTA& ha, EFTA& hb);

  inline EFTA& getEx() { return getEMFInstance()->getEx(); }
  inline EFTA& getEy() { return getEMFInstance()->getEy(); }
  inline EFTA& getEz() { return getEMFInstance()->getEz(); }
  inline EFTA& getHx() { return getEMFInstance()->getHx(); }
  inline EFTA& getHy() { return getEMFInstance()->getHy(); }
  inline EFTA& getHz() { return getEMFInstance()->getHz(); }
};
}  // namespace xfdtd
#endif  // _PERFECT_MATCH_LAYER_H_
