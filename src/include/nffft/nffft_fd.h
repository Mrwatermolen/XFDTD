#ifndef _NFFFT_FD_H_
#define _NFFFT_FD_H_

#include <cstddef>
#include <memory>
#include <utility>

#include "electromagnetic_field/electromagnetic_field.h"
#include "mesh/grid_box.h"
#include "nffft/nffft.h"
#include "util/type_define.h"

namespace xfdtd {
class NffftFd : public NFFFT {
 public:
  // define this constructor in header file
  NffftFd(SpatialIndex distance_x, SpatialIndex distance_y,
          SpatialIndex distance_z, xt::xarray<double> frequencies,
          xt::xarray<double> theta, xt::xarray<double> phi,
          std::filesystem::path output_dir_path)
      : NFFFT(distance_x, distance_y, distance_z, std::move(output_dir_path)),
        _frequencies{std::move(frequencies)},
        _wave_number{2 * constant::PI * _frequencies / constant::C_0},
        _theta{std::move(theta)},
        _phi{std::move(phi)},
        _number_frequecies{_frequencies.size()},
        _number_theta{_theta.size()},
        _number_phi{_phi.size()} {}
  NffftFd(const NffftFd &) = delete;
  NffftFd(NffftFd &&) = delete;
  NffftFd &operator=(const NffftFd &) = delete;
  NffftFd &operator=(NffftFd &&) = delete;
  ~NffftFd() override = default;

  void init(std::unique_ptr<GridBox> output_box, std::shared_ptr<EMF> emf,
            size_t total_time_steps, double dt, double dx, double dy,
            double dz) override;
  void update(size_t current_time_step) override;
  void outputData() override;

 private:
  xt::xarray<double> _frequencies, _wave_number, _theta, _phi;
  size_t _number_frequecies, _number_theta, _number_phi;

  EFFA _e_theta, _e_phi, _h_theta, _h_phi;
  EFFA _f_theta, _f_phi, _a_theta, _a_phi;

  EFFA _a_x, _a_y, _a_z;
  EFFA _f_x, _f_y, _f_z;

  EFFA _jx_yn, _jx_yp, _jx_zn, _jx_zp;
  EFFA _jy_xn, _jy_xp, _jy_zn, _jy_zp;
  EFFA _jz_xn, _jz_xp, _jz_yn, _jz_yp;
  EFFA _mx_yn, _mx_yp, _mx_zn, _mx_zp;
  EFFA _my_xn, _my_xp, _my_zn, _my_zp;
  EFFA _mz_xn, _mz_xp, _mz_yn, _mz_yp;

  EFFA _frequecy_transform_j;
  EFFA _frequecy_transform_m;

  void initDFT();

  void updateXN(size_t current_time_step);
  void updateXP(size_t current_time_step);
  void updateYN(size_t current_time_step);
  void updateYP(size_t current_time_step);
  void updateZN(size_t current_time_step);
  void updateZP(size_t current_time_step);

  void calculateFarfield();
  void calculateFarfieldX(int n, int t, int p, double sin_t_cos_p,
                          double sin_t_sin_p, double cos_t);
  void calculateFarfieldY(int n, int t, int p, double sin_t_cos_p,
                          double sin_t_sin_p, double cos_t);
  void calculateFarfieldZ(int n, int t, int p, double sin_t_cos_p,
                          double sin_t_sin_p, double cos_t);

  void calculateFarfield1(const PointVector &r, SpatialIndex range_a_start,
                          SpatialIndex range_a_end, SpatialIndex range_b_start,
                          SpatialIndex range_b_end, SpatialIndex range_c_start,
                          SpatialIndex range_c_end, const PointVector &offset);

  void outputFarFieldParamertes();
};
}  // namespace xfdtd

#endif  // _NFFFT_FD_H_