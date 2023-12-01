#ifndef _NFFFT_2D_TEST_H_
#define _NFFFT_2D_TEST_H_

#include "nffft/nffft.h"
#include "util/constant.h"
#include "util/type_define.h"
namespace xfdtd {
class NFFFT2DTEST : public NFFFT {
 public:
  NFFFT2DTEST(SpatialIndex distance_x, SpatialIndex distance_y,
              SpatialIndex direction_z, double far_theta, double far_phi,
              std::string output_dir_path);
  ~NFFFT2DTEST() override = default;

  void update() override;
  void outputData() override;

  void init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<const FDTDBasicCoff> fdtd_basic_coff,
            std::shared_ptr<const EMF> emf) override;

 private:
  EFTA _jx_yn, _jx_yp, _jx_zn, _jx_zp;
  EFTA _jy_xn, _jy_xp, _jy_zn, _jy_zp;
  EFTA _jz_xn, _jz_xp, _jz_yn, _jz_yp;
  EFTA _mx_yn, _mx_yp, _mx_zn, _mx_zp;
  EFTA _my_xn, _my_xp, _my_zn, _my_zp;
  EFTA _mz_xn, _mz_xp, _mz_yn, _mz_yp;
  double _farfield_phi{constant::PI * 1.25};
  double _farfield_theta{constant::PI / 2};

  void updateXN(size_t current_time_step);
  void updateXP(size_t current_time_step);
  void updateYN(size_t current_time_step);
  void updateYP(size_t current_time_step);
  void updateZN(size_t current_time_step);
  void updateZP(size_t current_time_step);

  void calculateFarfield();

  double getFarfieldPhi() { return _farfield_phi; }
};
}  // namespace xfdtd

#endif  // _NFFFT_2D_TEST_H_