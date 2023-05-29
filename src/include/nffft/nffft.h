#ifndef _NFFFT_H_
#define _NFFFT_H_

#include <memory>
#include <tuple>
#include <utility>

#include "electromagnetic_field/electromagnetic_field.h"
#include "mesh/grid_box.h"
#include "util/constant.h"
#include "util/type_define.h"

namespace xfdtd {
class NFFFT {
 public:
  NFFFT(SpatialIndex distance_x, SpatialIndex distance_y,
        SpatialIndex direction_z, double far_tehta, double far_phi);

  inline size_t getTotalTimeSteps() const { return _total_time_steps; }
  inline size_t getCurrentTimeStep() const { return _current_time_step; }
  inline std::tuple<SpatialIndex, SpatialIndex, SpatialIndex> getDistance() {
    return {_distance_x, _distance_y, _distance_z};
  }
  inline SpatialIndex getOutputBoundaryStartX() const {
    return _output_box->getStartIndexX();
  }
  inline SpatialIndex getOutputBoundaryEndX() const {
    return _output_box->getEndIndexX();
  }
  inline SpatialIndex getOutputBoundaryStartY() const {
    return _output_box->getStartIndexY();
  }
  inline SpatialIndex getOutputBoundaryEndY() const {
    return _output_box->getEndIndexY();
  }
  inline SpatialIndex getOutputBoundaryStartZ() const {
    return _output_box->getStartIndexZ();
  }
  inline SpatialIndex getOutputBoundaryEndZ() const {
    return _output_box->getEndIndexZ();
  }
  inline SpatialIndex getOutputBoundaryNx() const {
    return _output_box->getNx();
  }
  inline SpatialIndex getOutputBoundaryNy() const {
    return _output_box->getNy();
  }
  inline SpatialIndex getOutputBoundaryNz() const {
    return _output_box->getNz();
  }
  inline SpatialIndex getOutputBoundaryCenterX() const {
    return (getOutputBoundaryEndX() + getOutputBoundaryStartX()) / 2;
  }
  inline SpatialIndex getOutputBoundaryCenterY() const {
    return (getOutputBoundaryEndY() + getOutputBoundaryStartY()) / 2;
  }
  inline SpatialIndex getOutputBoundaryCenterZ() const {
    return (getOutputBoundaryEndZ() + getOutputBoundaryStartZ()) / 2;
  }
  inline const std::shared_ptr<const EMF> &getEMFInstance() const {
    return _emf;
  }
  inline double getFarfieldPhi() const { return _farfield_phi; }
  inline double getFarfieldTheta() const { return _farfield_theta; }
  inline double getDt() const { return _dt; }
  inline double getDx() const { return _dx; }
  inline double getDy() const { return _dy; }

  // WARNING: dz is always 1
  inline double getDz() const { return 1; }

  void update(size_t current_time_step);
  void outputData();

  void init(std::unique_ptr<GridBox> _output_box, std::shared_ptr<EMF> emf,
            size_t total_time_steps, double dt, double dx, double dy,
            double dz);

 private:
  SpatialIndex _distance_x, _distance_y, _distance_z;
  std::unique_ptr<GridBox> _output_box;
  std::shared_ptr<const EMF> _emf;
  size_t _total_time_steps;
  DoubleArrary3D _jmy_xn;
  DoubleArrary3D _jmy_xp;
  DoubleArrary3D _jmx_yn;
  DoubleArrary3D _jmx_yp;
  DoubleArrary3D _jez_xn;
  DoubleArrary3D _jez_xp;
  DoubleArrary3D _jez_yn;
  DoubleArrary3D _jez_yp;
  size_t _current_time_step{0};
  double _farfield_phi{constant::PI * 1.25};
  double _farfield_theta{constant::PI / 2};
  double _dt;
  double _dx;
  double _dy;
  double _dz;

  void getJmXN();
  void getJeJmXN();
  void getJeJmYN();
  void getJeJmXP();
  void getJeJmYP();

  void caculateFarfield();
};
}  // namespace xfdtd

#endif  // _NFFFT_H_