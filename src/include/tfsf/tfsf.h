#ifndef _TFSF_H_
#define _TFSF_H_

#include <memory>

#include "electromagnetic_field/electromagnetic_field.h"
#include "grid/grid_box.h"
#include "shape/cube.h"
#include "util/type_define.h"
#include "waveform/waveform.h"
namespace xfdtd {
/**
 * @brief The total field/scattered field (TFSF) boundary condition is a
 * technique for introducing far-zone sources (plane wave) into the FDTD
 * computational domain. Introduction of plane waves using the IFA(1D Incident
 * Field Array) method.
 *
 */
class TFSF {
 public:
  /**
   * @brief Construct a new TFSF object
   *
   * @param distance_x The distance in x between the TFSF boundary and the
   * computational domain boundary.
   * @param distance_y The distance in y between the TFSF boundary and the
   * computational domain boundary.
   * @param distance_z The distance in z between the TFSF boundary and the
   * computational domain boundary.
   * @param e_0 The amplitude of the incident field.
   * @param theta_inc
   * @param phi_inc
   * @param psi The incident wave polarization angle.
   * @param waveform The waveform of the incident field.
   */
  TFSF(SpatialIndex distance_x, SpatialIndex distance_y,
       SpatialIndex distance_z, double e_0, double theta_inc, double phi_inc,
       double psi, std::unique_ptr<Waveform> waveform);
  TFSF(TFSF&& ohter) = default;
  virtual ~TFSF() = default;

  inline std::tuple<SpatialIndex, SpatialIndex, SpatialIndex> getDistance()
      const {
    return {_distance_x, _distance_y, _distance_z};
  }
  inline double getIncidentFieldAmplitude() const { return _e_0; }
  inline double getIncidentTheta() const { return _theta_inc; }
  inline double getIncidentPhi() const { return _phi_inc; }
  inline double getIncidentPsi() const { return _psi; }

  inline double getIncidentFieldWaveformValueByTime(double time) {
    return _e_0 * _waveform->getValueByTime(time);
  }

  inline PointVector getKVector() const { return _k; }
  inline SpatialIndex getStartIndexX() const {
    return _tfsf_grid_box->getGridStartIndexX();
  }
  inline SpatialIndex getEndIndexX() const {
    return getStartIndexX() + getNx();
  }
  inline SpatialIndex getStartIndexY() const {
    return _tfsf_grid_box->getGridStartIndexY();
  }
  inline SpatialIndex getEndIndexY() const {
    return getStartIndexY() + getNy();
  }
  inline SpatialIndex getStartIndexZ() const {
    return _tfsf_grid_box->getGridStartIndexZ();
  }
  inline SpatialIndex getEndIndexZ() const {
    return getStartIndexZ() + getNz();
  }
  inline SpatialIndex getNx() const { return _tfsf_grid_box->getGridNumX(); }
  inline SpatialIndex getNy() const { return _tfsf_grid_box->getGridNumY(); }
  inline SpatialIndex getNz() const { return _tfsf_grid_box->getGridNumZ(); }

  inline PointVector getKInc() const { return _k; }

  void setEMFInstance(std::shared_ptr<EMF> emf) { _emf = std::move(emf); };

  virtual void init(double dx, double dy, double dz, double dt,
                    std::unique_ptr<GridBox> tfsf_grid_box) = 0;
  virtual void updateIncidentField(size_t current_time_step) = 0;
  virtual void updateH() = 0;
  /**
   * @brief correct E feild in TFSF boundary
   *
   */
  virtual void updateE() = 0;

 protected:
  void defaultInitTFSF(double dx, double dy, double dz, double dt,
                       std::unique_ptr<GridBox> tfsf_grid_box);

  inline double getIncidentSinTheta() const { return _sin_theta_inc; }
  inline double getIncidentCosTheta() const { return _cos_theta_inc; }
  inline double getIncidentSinPhi() const { return _sin_phi_inc; }
  inline double getIncidentCosPhi() const { return _cos_phi_inc; }
  inline double getIncidentSinPsi() const { return _sin_psi; }
  inline double getIncidentCosPsi() const { return _cos_psi; }
  inline double getDt() const { return _dt; }
  inline double getDx() const { return _dx; }
  inline double getDy() const { return _dy; }
  inline double getDz() const { return _dz; }

  inline double& getEx(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getEx(i, j, k);
  }
  inline double& getEy(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getEy(i, j, k);
  }
  inline double& getEz(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getEz(i, j, k);
  }
  inline double& getHx(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getHx(i, j, k);
  }
  inline double& getHy(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getHy(i, j, k);
  }
  inline double& getHz(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getHz(i, j, k);
  }

 private:
  SpatialIndex _distance_x;
  SpatialIndex _distance_y;
  SpatialIndex _distance_z;
  double _e_0;
  double _theta_inc;
  double _phi_inc;
  double _psi;
  double _sin_theta_inc;
  double _cos_theta_inc;
  double _sin_phi_inc;
  double _cos_phi_inc;
  double _sin_psi;
  double _cos_psi;
  PointVector _k;
  std::shared_ptr<Waveform> _waveform;

  double _dt;
  double _dx;
  double _dy;
  double _dz;
  std::unique_ptr<GridBox> _tfsf_grid_box;
  std::shared_ptr<EMF> _emf;
};
}  // namespace xfdtd

#endif  // _TFSF_H_
