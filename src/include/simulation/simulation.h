#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <cstddef>
#include <filesystem>
#include <memory>
#include <vector>

#include "boundary/boundary.h"
#include "electromagnetic_field/electromagnetic_field.h"
#include "fdtd_basic_coff/fdtd_basic_coff.h"
#include "grid/grid_space.h"
#include "lumped_element/lumped_element.h"
#include "monitor/monitor.h"
#include "network/network.h"
#include "nffft/nffft.h"
#include "object/object.h"
#include "tfsf/tfsf.h"
#include "util/type_define.h"

namespace xfdtd {

using ObjectArray = std::vector<std::shared_ptr<Object>>;
using BoundaryArray = std::vector<std::shared_ptr<Boundary>>;
using MonitorArray = std::vector<std::shared_ptr<Monitor>>;

class Simulation {
 public:
  inline static constexpr float DEFAULT_CFL{0.99};

  explicit Simulation(double cell_size, float cfl = DEFAULT_CFL);
  Simulation(double cell_size, ObjectArray objects, BoundaryArray boundaries,
             std::unique_ptr<TFSF> tfsf, float cfl = DEFAULT_CFL);
  Simulation(double cell_size, ObjectArray objects, BoundaryArray boundaries,
             std::unique_ptr<TFSF> tfsf, std::unique_ptr<NFFFT> nffft,
             float cfl = DEFAULT_CFL);
  Simulation(double cell_size, ObjectArray objects,
             std::vector<std::shared_ptr<LumpedElement>> lumped_elements,
             BoundaryArray boundaries = {}, MonitorArray monitor = {},
             float cfl = DEFAULT_CFL);

  Simulation(double cell_size, ObjectArray objects,
             std::vector<std::shared_ptr<LumpedElement>> lumped_elements,
             MonitorArray monitors = {}, float cfl = DEFAULT_CFL);

  Simulation(double cell_size, ObjectArray objects, MonitorArray monitors = {},
             float cfl = DEFAULT_CFL);

  Simulation(double cell_size, ObjectArray objects,
             std::vector<std::shared_ptr<LumpedElement>> lumped_elements,
             BoundaryArray boundaries, MonitorArray monitors,
             std::unique_ptr<Network> network, float cfl = DEFAULT_CFL);

  void addObject(std::unique_ptr<Object> object);
  void addTFSFSource(std::unique_ptr<TFSF> tfsf);
  void addNFFFT(std::unique_ptr<NFFFT> nffft);
  void addMonitor(std::unique_ptr<Monitor> monitor);

  void checkRun(size_t total_time_steps);
  void run(size_t total_time_steps);

  // inline double getDx() const { return _grid_space->getDx(); }
  // inline double getDy() const { return _grid_space->getDy(); }
  // inline double getDz() const { return _grid_space->getDz(); }
  inline double getDt() const { return _fdtd_basic_coff->getDt(); }
  inline SpatialIndex getNx() const { return _grid_space->getGridNumX(); }
  inline SpatialIndex getNy() const { return _grid_space->getGridNumY(); }
  inline SpatialIndex getNz() const { return _grid_space->getGridNumZ(); }
  inline EFTA& getCexe() { return _fdtd_basic_coff->getCexe(); }
  inline EFTA& getCexhy() { return _fdtd_basic_coff->getCexhy(); }
  inline EFTA& getCexhz() { return _fdtd_basic_coff->getCexhz(); }
  inline EFTA& getCeye() { return _fdtd_basic_coff->getCeye(); }
  inline EFTA& getCeyhz() { return _fdtd_basic_coff->getCeyhz(); }
  inline EFTA& getCeyhx() { return _fdtd_basic_coff->getCeyhx(); }
  inline EFTA& getCeze() { return _fdtd_basic_coff->getCeze(); }
  inline EFTA& getCezhx() { return _fdtd_basic_coff->getCezhx(); }
  inline EFTA& getCezhy() { return _fdtd_basic_coff->getCezhy(); }
  inline EFTA& getChxh() { return _fdtd_basic_coff->getChxh(); }
  inline EFTA& getChxey() { return _fdtd_basic_coff->getChxey(); }
  inline EFTA& getChxez() { return _fdtd_basic_coff->getChxez(); }
  inline EFTA& getChyh() { return _fdtd_basic_coff->getChyh(); }
  inline EFTA& getChyez() { return _fdtd_basic_coff->getChyez(); }
  inline EFTA& getChyex() { return _fdtd_basic_coff->getChyex(); }
  inline EFTA& getChzh() { return _fdtd_basic_coff->getChzh(); }
  inline EFTA& getChzex() { return _fdtd_basic_coff->getChzex(); }
  inline EFTA& getChzey() { return _fdtd_basic_coff->getChzey(); }

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

  void outputTFSFIncidentWaveFastFourierTransform(const std::filesystem::path &path);

 private:
  ObjectArray _objects;
  BoundaryArray _boundaries;
  std::vector<std::shared_ptr<LumpedElement>> _lumped_elements;
  std::unique_ptr<TFSF> _tfsf;
  std::unique_ptr<NFFFT> _nffft;
  MonitorArray _monitors;
  std::unique_ptr<Network> _network;
  std::shared_ptr<GridSpace> _grid_space;
  std::shared_ptr<EMF> _emf;
  std::shared_ptr<FDTDBasicCoff> _fdtd_basic_coff;

  size_t _total_time_steps;
  size_t _current_time_step;
  size_t _nx, _ny, _nz;

  bool _is_exist_dispersive_material{false};
  xt::xarray<bool> _is_exist_dispersive_material_array;

  void init();
  // void initMaterialGrid();
  void initGridSpace();
  void initFDTDBasicCoff();
  void initEMInstance();
  void initObject();
  void initTFSF();
  void initNFFFT();
  void initUpdateCoefficient();
  void initBoundaryCondition();
  void initMonitor();
  void initLumpedElement();
  void initNetwork();

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

  void correctE();

  void correctH();

  void outputData();

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
