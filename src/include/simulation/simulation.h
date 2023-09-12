#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <cstddef>
#include <memory>
#include <vector>

#include "boundary/boundary.h"
#include "electromagnetic_field/electromagnetic_field.h"
#include "monitor/monitor.h"
#include "nffft/nffft.h"
#include "object/object.h"
#include "shape/cube.h"
#include "tfsf/tfsf.h"
#include "util/type_define.h"

namespace xfdtd {
class TFSF2D;
class Simulation {
 public:
  inline static constexpr float DEFAULT_CFL{0.99};
  Simulation() : _emf{std::make_shared<EMF>()} {}

  explicit Simulation(double cell_size, float cfl = DEFAULT_CFL);
  Simulation(double cell_size, ObjectArray objects, BoundaryArray boundaries,
             std::unique_ptr<TFSF> tfsf, float cfl = DEFAULT_CFL);
  Simulation(double cell_size, ObjectArray objects, BoundaryArray boundaries,
             std::unique_ptr<TFSF> tfsf, std::unique_ptr<NFFFT> nffft,
             float cfl = DEFAULT_CFL);

  void setCellSize(double cell_size);
  void addObject(std::unique_ptr<Object> object);
  void addTFSFSource(std::unique_ptr<TFSF> tfsf);
  void addNFFFT(std::unique_ptr<NFFFT> nffft);
  void addMonitor(std::unique_ptr<Monitor> monitor);

  void checkRun(size_t time_steps);
  void run(size_t time_steps);

  inline double getDx() const { return _dx; }
  inline double getDy() const { return _dy; }
  inline double getDz() const { return _dz; }
  inline double getDt() const { return _dt; }
  inline SpatialIndex getNx() const { return _nx; }
  inline SpatialIndex getNy() const { return _ny; }
  inline SpatialIndex getNz() const { return _nz; }
  inline EFTA& getCexe() { return _cexe; }
  inline EFTA& getCexhy() { return _cexhy; }
  inline EFTA& getCexhz() { return _cexhz; }
  inline EFTA& getCeye() { return _ceye; }
  inline EFTA& getCeyhz() { return _ceyhz; }
  inline EFTA& getCeyhx() { return _ceyhx; }
  inline EFTA& getCeze() { return _ceze; }
  inline EFTA& getCezhx() { return _cezhx; }
  inline EFTA& getCezhy() { return _cezhy; }
  inline EFTA& getChxh() { return _chxh; }
  inline EFTA& getChxey() { return _chxey; }
  inline EFTA& getChxez() { return _chxez; }
  inline EFTA& getChyh() { return _chyh; }
  inline EFTA& getChyez() { return _chyez; }
  inline EFTA& getChyex() { return _chyex; }
  inline EFTA& getChzh() { return _chzh; }
  inline EFTA& getChzex() { return _chzex; }
  inline EFTA& getChzey() { return _chzey; }

  inline std::shared_ptr<EMF> getEMFInstance() { return _emf; }
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

 private:
  // simulation parameter
  double _dx;
  double _dy;
  double _dz;
  float _cfl{DEFAULT_CFL};

  ObjectArray _objects;
  BoundaryArray _boundaries;
  std::unique_ptr<TFSF> _tfsf;
  std::unique_ptr<NFFFT> _nffft;
  MonitorArray _monitors;

  SpatialIndex _nx;
  SpatialIndex _ny;
  SpatialIndex _nz;
  double _dt;
  size_t _time_steps;
  size_t _current_time_step;

  xt::xarray<std::shared_ptr<YeeCell>> _grid_space;
  std::unique_ptr<Cube> _simulation_box;

  std::shared_ptr<EMF> _emf;

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

  bool _is_exist_dispersive_material{false};
  xt::xarray<bool> _is_exist_dispersive_material_array;

  void init();
  void initMaterialGrid();
  void initTFSF();
  void initNFFFT();
  void initUpdateCoefficient();
  void initBoundaryCondition();
  // TODO(franzero): fix it
  void initMonitor();

  void calculateDomainSize();
  void gridSimulationSpace();
  void allocateArray();

  void updateTFSFIncidentField();
  void updateH();
  void updateTFSFH();
  void updateBoundaryH();
  void updateE();
  void updateTFSFE();
  void updateBoundaryE();
  void updateNFFFT();
  void updateMonitor();
  void updateEWithDispersiveMaterial();

  inline void allocateEx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz) {
    _emf->allocateEx(nx, ny, nz);
  }

  inline void allocateEy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz) {
    _emf->allocateEy(nx, ny, nz);
  }

  inline void allocateEz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz) {
    _emf->allocateEz(nx, ny, nz);
  }

  inline void allocateHx(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz) {
    _emf->allocateHx(nx, ny, nz);
  }

  inline void allocateHy(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz) {
    _emf->allocateHy(nx, ny, nz);
  }

  inline void allocateHz(SpatialIndex nx, SpatialIndex ny, SpatialIndex nz) {
    _emf->allocateHz(nx, ny, nz);
  }
};

}  // namespace xfdtd

#endif  // _SIMULATION_H_
