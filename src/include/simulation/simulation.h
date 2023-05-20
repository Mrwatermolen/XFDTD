#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <cstddef>
#include <memory>
#include <vector>

#include "additive_source/source.h"
#include "boundary/boundary.h"
#include "electromagnetic.h"
#include "monitor/monitor.h"
#include "object/object.h"
#include "shape/cube.h"
#include "simulation/yee_cell.h"
#include "tfsf/tfsf.h"
#include "util/type_define.h"

namespace xfdtd {

class Simulation {
 public:
  inline static constexpr float DEFAULT_CFL{0.99};
  Simulation(double cell_size, ObjectArray objects, SourceArray sources,
             std::unique_ptr<TFSF> tfsf, BoundaryArray boundaries = {},
             MonitorArray _monitors = {}, float cfl = DEFAULT_CFL);

  void checkRun(size_t time_steps);
  void run(size_t time_steps);
  void output();

 private:
  // simulation parameter
  double _dx;
  double _dy;
  double _dz;
  float _cfl{DEFAULT_CFL};

  ObjectArray _objects;
  SourceArray _sources;
  BoundaryArray _boundaries;
  std::unique_ptr<TFSF> _tfsf;
  MonitorArray _monitors;

  SpatialIndex _nx;
  SpatialIndex _ny;
  SpatialIndex _nz;
  double _dt;
  size_t _time_steps;
  size_t _current_time_step;
  std::vector<double> _time_array;  // doubt that it is necessary

  YeeCellArray _grid_space;
  std::unique_ptr<Cube> _simulation_box;

  EFTA _cexe;
  EFTA _cexhy;
  EFTA _cexhz;
  EFTA _cexje;
  EFTA _ceye;
  EFTA _ceyhz;
  EFTA _ceyhx;
  EFTA _ceyje;
  EFTA _ceze;
  EFTA _cezhx;
  EFTA _cezhy;
  EFTA _cezje;

  EFTA _chxh;
  EFTA _chxey;
  EFTA _chxez;
  EFTA _chxjm;
  EFTA _chyh;
  EFTA _chyez;
  EFTA _chyex;
  EFTA _chyjm;
  EFTA _chzh;
  EFTA _chzex;
  EFTA _chzey;
  EFTA _chzjm;

  EFTA _eps_x;
  EFTA _eps_y;
  EFTA _eps_z;
  EFTA _sigma_e_x;
  EFTA _sigma_e_y;
  EFTA _sigma_e_z;

  EFTA _mu_x;
  EFTA _mu_y;
  EFTA _mu_z;
  EFTA _sigma_m_x;
  EFTA _sigma_m_y;
  EFTA _sigma_m_z;

  void init();
  void initMaterialGrid();
  void initSource();
  void initTFSF();
  void initUpdateCoefficient();
  void initBondaryCondition();
  void initMonitor();

  void caculateDomainSize();
  void gridSimualtionSpace();
  void allocateArray();
  void caculateMaterialComponent();

  void updateAddtiveSource();
  void updateTFSFIncidentField();
  void updateH();
  void updateTFSFH();
  void updateBoundaryH();
  void updateE();
  void updateTFSFE();
  void updateBoundaryE();
  void updateMonitor();

  void handleHardPointSource(Source* source);
  void handlePMLBoundaryH(std::shared_ptr<Boundary>& boundary);
  void handlePMLBoundaryE(std::shared_ptr<Boundary>& boundary);
};

}  // namespace xfdtd

#endif  // _SIMULATION_H_
