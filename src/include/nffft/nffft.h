#ifndef _NFFFT_H_
#define _NFFFT_H_

#include <filesystem>
#include <memory>
#include <tuple>

#include "electromagnetic_field/electromagnetic_field.h"
#include "fdtd_basic_coff/fdtd_basic_coff.h"
#include "grid/grid_box.h"
#include "grid/grid_space.h"
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

  virtual void update() = 0;
  virtual void outputData() = 0;
  virtual void init(std::shared_ptr<const GridSpace> grid_space,
                    std::shared_ptr<const FDTDBasicCoff> fdtd_basic_coff,
                    std::shared_ptr<const EMF> emf) = 0;

 protected:
  void defaultInit(std::shared_ptr<const GridSpace> grid_space,
                   std::shared_ptr<const FDTDBasicCoff> fdtd_basic_coff,
                   std::shared_ptr<const EMF> emf);

  inline size_t getTotalTimeSteps() {
    return _fdtd_basic_coff->getTotalTimeStep();
  }
  size_t getCurrentTimeStep() { return _fdtd_basic_coff->getCurrentTimeStep(); }
  inline double getDt() const { return _fdtd_basic_coff->getDt(); }
  inline double getDx() const { return _grid_space->getGridBaseSizeX(); }
  inline double getDy() const { return _grid_space->getGridBaseSizeY(); }
  inline double getDz() const { return _grid_space->getGridBaseSizeZ(); }

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

 protected:
  const GridSpace *getGridSpace() const { return _grid_space.get(); }

  const FDTDBasicCoff *getFDTDBasicCoff() const {
    return _fdtd_basic_coff.get();
  }

  const EMF *getEMFInstance() const { return _emf.get(); }

 private:
  SpatialIndex _distance_x, _distance_y, _distance_z;
  std::filesystem::path _output_dir_path;

  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const FDTDBasicCoff> _fdtd_basic_coff;
  std::shared_ptr<const EMF> _emf;
  std::unique_ptr<GridBox> _output_box;
};
}  // namespace xfdtd

#endif  // _NFFFT_H_