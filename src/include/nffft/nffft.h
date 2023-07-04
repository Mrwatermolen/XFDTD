#ifndef _NFFFT_H_
#define _NFFFT_H_

#include <filesystem>
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
        SpatialIndex direction_z, std::filesystem::path output_dir_path);
  virtual ~NFFFT() = default;

  inline std::tuple<SpatialIndex, SpatialIndex, SpatialIndex> getDistance() {
    return {_distance_x, _distance_y, _distance_z};
  }

  inline SpatialIndex getStartIndexX() { return _output_box->getStartIndexX(); }
  inline SpatialIndex getStartIndexY() { return _output_box->getStartIndexY(); }
  inline SpatialIndex getStartIndexZ() { return _output_box->getStartIndexZ(); }
  inline SpatialIndex getEndIndexX() { return _output_box->getEndIndexX(); }
  inline SpatialIndex getEndIndexY() { return _output_box->getEndIndexY(); }
  inline SpatialIndex getEndIndexZ() { return _output_box->getEndIndexZ(); }
  inline SpatialIndex getNx() { return _output_box->getNx(); }
  inline SpatialIndex getNy() { return _output_box->getNy(); }
  inline SpatialIndex getNz() { return _output_box->getNz(); }
  inline size_t getTotalTimeSteps() { return _total_time_steps; }
  inline double getDt() const { return _dt; }
  inline double getDx() const { return _dx; }
  inline double getDy() const { return _dy; }
  inline double getDz() const { return _dz; }
  inline std::filesystem::path getOutputDirPath() { return _output_dir_path; }

  virtual void update(size_t current_time_step) = 0;
  virtual void outputData() = 0;

  virtual void init(std::unique_ptr<GridBox> output_box,
                    std::shared_ptr<EMF> emf, size_t total_time_steps,
                    double dt, double dx, double dy, double dz) = 0;

 protected:
  void defaultInit(std::unique_ptr<GridBox> output_box,
                   std::shared_ptr<EMF> emf, size_t total_time_steps, double dt,
                   double dx, double dy, double dz);

  double getEx(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getEx(i, j, k);
  }
  double getEy(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getEy(i, j, k);
  }
  double getEz(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getEz(i, j, k);
  }
  double getHx(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getHx(i, j, k);
  }
  double getHy(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getHy(i, j, k);
  }
  double getHz(SpatialIndex i, SpatialIndex j, SpatialIndex k) {
    return _emf->getHz(i, j, k);
  }

 private:
  SpatialIndex _distance_x, _distance_y, _distance_z;
  std::filesystem::path _output_dir_path;

  std::unique_ptr<GridBox> _output_box;
  std::shared_ptr<const EMF> _emf;
  size_t _total_time_steps;
  double _dt;
  double _dx;
  double _dy;
  double _dz;
};
}  // namespace xfdtd

#endif  // _NFFFT_H_