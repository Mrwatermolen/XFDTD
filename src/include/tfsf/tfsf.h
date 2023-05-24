#ifndef _TFSF_H_
#define _TFSF_H_

#include <memory>

#include "electromagnetic_field/electromagnetic_field.h"
#include "shape/cube.h"
#include "util/type_define.h"
#include "waveform/waveform.h"
namespace xfdtd {
class TFSF {
 public:
  TFSF(SpatialIndex distance_x, SpatialIndex distance_y,
       SpatialIndex distance_z, double theta_inc, double phi_inc,
       double e_theta, double e_phi, std::unique_ptr<Waveform> waveform);
  TFSF(TFSF&& ohter) = default;
  virtual ~TFSF() = default;

  struct TFSFBoundaryIndex {
    SpatialIndex start_x;
    SpatialIndex start_y;
    SpatialIndex start_z;
    SpatialIndex nx;
    SpatialIndex ny;
    SpatialIndex nz;
  };

  inline std::tuple<SpatialIndex, SpatialIndex, SpatialIndex> getDistance()
      const {
    return {_distance_x, _distance_y, _distance_z};
  }
  inline double getIncidentTheta() const { return _theta_inc; }
  inline double getIncidentPhi() const { return _phi_inc; }
  inline double getETehta() const { return _e_theta; }
  inline double getEPhi() const { return _e_phi; }
  inline double getIncidentFieldWaveformValueByTime(double time) {
    return _waveform->getValueByTime(time);
  }
  inline double getDt() const { return _dt; }
  inline double getDx() const { return _dx; }
  inline double getDy() const { return _dy; }
  inline double getDz() const { return _dz; }

  /**
   * @brief get the wave vector of the incident field
   */
  inline PointVector getKVector() const { return _k; }

  /**
   * @brief Get the index of the TFSF boundary with normal vector xn
   */
  inline SpatialIndex getStartIndexX() const {
    return _tfsf_boundary_index.start_x;
  }
  /**
   * @brief Get the index of the TFSF boundary with normal vector xp
   */
  inline SpatialIndex getEndIndexX() const {
    return getStartIndexX() + getNx();
  }
  /**
   * @brief Get the index of the TFSF boundary with normal vector yn
   */
  inline SpatialIndex getStartIndexY() const {
    return _tfsf_boundary_index.start_y;
  }
  /**
   * @brief Get the index of the TFSF boundary with normal vector yp
   */
  inline SpatialIndex getEndIndexY() const {
    return getStartIndexX() + getNx();
  }
  inline SpatialIndex getStartIndexZ() const {
    return _tfsf_boundary_index.start_z;
  }
  inline SpatialIndex getEndIndexZ() const {
    return getStartIndexZ() + getNz();
  }
  /**
   * @brief Get the number of cells in the x direction of the TFSF boundary
   */
  inline SpatialIndex getNx() const { return _tfsf_boundary_index.nx; }
  /**
   * @brief Get the number of cells in the y direction of the TFSF boundary
   */
  inline SpatialIndex getNy() const { return _tfsf_boundary_index.ny; }
  /**
   * @brief Get the number of cells in the z direction of the TFSF boundary
   */
  inline SpatialIndex getNz() const { return _tfsf_boundary_index.nz; }

  inline const Cube* getTFSFCubeBox() const {
    if (_tfsf_box == nullptr) {
      throw std::runtime_error("TFSF box is not initialized");
    }
    return _tfsf_box.get();
  }

  inline EFTA& getEx() { return _emf->getEx(); }
  inline EFTA& getEy() { return _emf->getEy(); }
  inline EFTA& getEz() { return _emf->getEz(); }
  inline EFTA& getHx() { return _emf->getHx(); }
  inline EFTA& getHy() { return _emf->getHy(); }
  inline EFTA& getHz() { return _emf->getHz(); }
  inline double& getEx(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getEx(i, j, k);
  }
  inline double& getExy(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
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

  void setEMFInstance(std::shared_ptr<EMF> emf) { _emf = std::move(emf); };

  virtual void init(const Cube* simulation_box, double dx, double dy, double dz,
                    double dt, TFSFBoundaryIndex tfsf_boundary_index) = 0;
  virtual void updateIncidentField(size_t current_time_step) = 0;
  virtual void updateH() = 0;
  virtual void updateE() = 0;

 protected:
  void defaultInitTFSF(const Cube* simulation_box, double dx, double dy,
                       double dz, double dt,
                       TFSFBoundaryIndex tfsf_boundary_index);

 private:
  SpatialIndex _distance_x;
  SpatialIndex _distance_y;
  SpatialIndex _distance_z;
  double _theta_inc;
  double _phi_inc;
  double _e_theta;
  double _e_phi;
  PointVector _k;
  std::unique_ptr<Waveform> _waveform;

  double _dt;
  double _dx;
  double _dy;
  double _dz;
  TFSFBoundaryIndex _tfsf_boundary_index;
  std::unique_ptr<Cube> _tfsf_box;
  std::shared_ptr<EMF> _emf;
};
}  // namespace xfdtd

#endif  // _TFSF_H_