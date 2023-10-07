#ifndef _XFDTD_FDTD_BASIC_COFF_H_
#define _XFDTD_FDTD_BASIC_COFF_H_

#include <cstddef>
#include <xtensor/xarray.hpp>

namespace xfdtd {

class FDTDBasicCoff {
 public:
  FDTDBasicCoff(double dx, double dy, double dz, double cfl);
  FDTDBasicCoff(double dx, double dy, double dz, double cfl,
                size_t total_time_step);

  void init(size_t nx, size_t ny, size_t nz, size_t total_time_step);

  void initCoff();

  double getDt() const;

  double getDx() const;

  double getDy() const;

  double getDz() const;

  size_t getTotalTimeStep() const;

  size_t getCurrentTimeStep() const;

  const xt::xarray<double>& getEpsX() const;

  const xt::xarray<double>& getEpsY() const;

  const xt::xarray<double>& getEpsZ() const;

  const xt::xarray<double>& getMuX() const;

  const xt::xarray<double>& getMuY() const;

  const xt::xarray<double>& getMuZ() const;

  const xt::xarray<double>& getSigmaX() const;

  const xt::xarray<double>& getSigmaY() const;

  const xt::xarray<double>& getSigmaZ() const;

  const xt::xarray<double>& getSigmaMx() const;

  const xt::xarray<double>& getSigmaMy() const;

  const xt::xarray<double>& getSigmaMz() const;

  const xt::xarray<double>& getCexe() const;

  const xt::xarray<double>& getCexhy() const;

  const xt::xarray<double>& getCexhz() const;

  const xt::xarray<double>& getCeye() const;

  const xt::xarray<double>& getCeyhz() const;

  const xt::xarray<double>& getCeyhx() const;

  const xt::xarray<double>& getCeze() const;

  const xt::xarray<double>& getCezhx() const;

  const xt::xarray<double>& getCezhy() const;

  const xt::xarray<double>& getChxh() const;

  const xt::xarray<double>& getChxey() const;

  const xt::xarray<double>& getChxez() const;

  const xt::xarray<double>& getChyh() const;

  const xt::xarray<double>& getChyez() const;

  const xt::xarray<double>& getChyex() const;

  const xt::xarray<double>& getChzh() const;

  const xt::xarray<double>& getChzex() const;

  const xt::xarray<double>& getChzey() const;

  xt::xarray<double>& getEpsX();

  xt::xarray<double>& getEpsY();

  xt::xarray<double>& getEpsZ();

  xt::xarray<double>& getMuX();

  xt::xarray<double>& getMuY();

  xt::xarray<double>& getMuZ();

  xt::xarray<double>& getSigmaX();

  xt::xarray<double>& getSigmaY();

  xt::xarray<double>& getSigmaZ();

  xt::xarray<double>& getSigmaMx();

  xt::xarray<double>& getSigmaMy();

  xt::xarray<double>& getSigmaMz();

  xt::xarray<double>& getCexe();

  xt::xarray<double>& getCexhy();

  xt::xarray<double>& getCexhz();

  xt::xarray<double>& getCeye();

  xt::xarray<double>& getCeyhz();

  xt::xarray<double>& getCeyhx();

  xt::xarray<double>& getCeze();

  xt::xarray<double>& getCezhx();

  xt::xarray<double>& getCezhy();

  xt::xarray<double>& getChxh();

  xt::xarray<double>& getChxey();

  xt::xarray<double>& getChxez();

  xt::xarray<double>& getChyh();

  xt::xarray<double>& getChyez();

  xt::xarray<double>& getChyex();

  xt::xarray<double>& getChzh();

  xt::xarray<double>& getChzex();

  xt::xarray<double>& getChzey();

  void setCurrentTimeStep(size_t current_time_step);

 private:
  double _cfl;
  double _dx, _dy, _dz;
  double _dt;
  size_t _total_time_step, _current_time_step;

  xt::xarray<double> _eps_x;
  xt::xarray<double> _eps_y;
  xt::xarray<double> _eps_z;
  xt::xarray<double> _mu_x;
  xt::xarray<double> _mu_y;
  xt::xarray<double> _mu_z;
  xt::xarray<double> _sigma_x;
  xt::xarray<double> _sigma_y;
  xt::xarray<double> _sigma_z;
  xt::xarray<double> _sigma_m_x;
  xt::xarray<double> _sigma_m_y;
  xt::xarray<double> _sigma_m_z;

  xt::xarray<double> _cexe;
  xt::xarray<double> _cexhy;
  xt::xarray<double> _cexhz;
  xt::xarray<double> _ceye;
  xt::xarray<double> _ceyhz;
  xt::xarray<double> _ceyhx;
  xt::xarray<double> _ceze;
  xt::xarray<double> _cezhx;
  xt::xarray<double> _cezhy;

  xt::xarray<double> _chxh;
  xt::xarray<double> _chxey;
  xt::xarray<double> _chxez;
  xt::xarray<double> _chyh;
  xt::xarray<double> _chyez;
  xt::xarray<double> _chyex;
  xt::xarray<double> _chzh;
  xt::xarray<double> _chzex;
  xt::xarray<double> _chzey;

  void allocateArray(size_t nx, size_t ny, size_t nz);
};

}  // namespace xfdtd

#endif  // _XFDTD_FDTD_BASIC_COFF_H_
