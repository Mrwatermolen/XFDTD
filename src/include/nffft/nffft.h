#ifndef _NFFFT_H_
#define _NFFFT_H_

#include <filesystem>
#include <memory>
#include <tuple>

#include "electromagnetic_field/electromagnetic_field.h"
#include "grid/grid_box.h"
#include "util/constant.h"
#include "util/type_define.h"

namespace xfdtd {
/**
 * @brief The near-field to far-field transformation (NFFFT) is a technique for
 * computing the far-field. The NFFFT is based on the discrete Fourier transform
 * (DFT) and the fast Fourier transform (FFT).
 *
 */
class NFFFT {
 public:
  /**
   * @brief Construct a new NFFFT object
   *
   * @param distance_x The distance in x between the output boundary and the
   * computational domain boundary
   * @param distance_y The distance in y between the output boundary and the
   * computational domain boundary
   * @param distance_z The distance in z between the output boundary and the
   * computational domain boundary
   * @param output_dir_path The path of the output directory
   */
  NFFFT(SpatialIndex distance_x, SpatialIndex distance_y,
        SpatialIndex distance_z, std::filesystem::path output_dir_path);
  NFFFT(const NFFFT &) = delete;
  NFFFT(NFFFT &&) = default;
  NFFFT &operator=(const NFFFT &) = delete;
  NFFFT &operator=(NFFFT &&) = default;
  virtual ~NFFFT() = default;

  inline std::tuple<SpatialIndex, SpatialIndex, SpatialIndex> getDistance() {
    return {_distance_x, _distance_y, _distance_z};
  }

  inline SpatialIndex getStartIndexX() {
    return _output_box->getGridStartIndexX();
  }
  inline SpatialIndex getStartIndexY() {
    return _output_box->getGridStartIndexY();
  }
  inline SpatialIndex getStartIndexZ() {
    return _output_box->getGridStartIndexZ();
  }
  inline SpatialIndex getEndIndexX() { return _output_box->getGridEndIndexX(); }
  inline SpatialIndex getEndIndexY() { return _output_box->getGridEndIndexY(); }
  inline SpatialIndex getEndIndexZ() { return _output_box->getGridEndIndexZ(); }
  inline SpatialIndex getNx() { return _output_box->getGridNumX(); }
  inline SpatialIndex getNy() { return _output_box->getGridNumY(); }
  inline SpatialIndex getNz() { return _output_box->getGridNumZ(); }
  inline SpatialIndex getCenterIndexX() {
    return _output_box->getGridCenterIndexX();
  }
  inline SpatialIndex getCenterIndexY() {
    return _output_box->getGridCenterIndexY();
  }
  inline SpatialIndex getCenterIndexZ() {
    return _output_box->getGridCenterIndexZ();
  }
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

  inline size_t getTotalTimeSteps() { return _total_time_steps; }
  inline double getDt() const { return _dt; }
  inline double getDx() const { return _dx; }
  inline double getDy() const { return _dy; }
  inline double getDz() const { return _dz; }
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