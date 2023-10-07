#ifndef _XFDTD_MONITOR_H_
#define _XFDTD_MONITOR_H_

#include <filesystem>
#include <memory>

#include "electromagnetic_field/electromagnetic_field.h"
#include "fdtd_basic_coff/fdtd_basic_coff.h"
#include "grid/grid_space.h"
#include "shape/shape.h"
namespace xfdtd {

class Monitor {
 public:
  Monitor() = default;

  Monitor(std::unique_ptr<Shape> shape, std::filesystem::path output_dir_path,
          std::string output_file_name);

  Monitor(const Monitor& other);

  Monitor& operator=(const Monitor& other);

  Monitor(Monitor&& other) noexcept = default;

  Monitor& operator=(Monitor&& other) noexcept = default;

  virtual ~Monitor() = default;

  virtual std::unique_ptr<Monitor> clone() const = 0;

  virtual void init(const std::shared_ptr<const FDTDBasicCoff>& fdtd_basic_coff,
                    const std::shared_ptr<const GridSpace>& grid_space,
                    const std::shared_ptr<const EMF>& emf) = 0;

  virtual void update() = 0;

  virtual void outputData() = 0;

  virtual const std::filesystem::path& getOutputPath() const;

  virtual const std::string& getOutputFileName() const;

  virtual void setOutputDirPath(const std::string& output_dir_path);

  virtual void setOutputFileName(const std::string& output_file_name);

 protected:
  void defaultInit(const std::shared_ptr<const FDTDBasicCoff>& fdtd_basic_coff,
                   const std::shared_ptr<const GridSpace>& grid_space,
                   const std::shared_ptr<const EMF>& emf);

  virtual const Shape* getShape() const;

  const GridSpace* getGridSpaceInstance() const;

  const FDTDBasicCoff* getFDTDBasicCoffInstance() const;

  const EMF* getEMFInstance() const;

 private:
  std::unique_ptr<Shape> _shape;
  std::filesystem::path _output_dir_path;
  std::string _output_file_name;
  std::shared_ptr<const FDTDBasicCoff> _fdtd_basic_coff;
  std::shared_ptr<const GridSpace> _grid_space;
  std::shared_ptr<const EMF> _emf;
};

}  // namespace xfdtd

#endif  // _XFDTD_MONITOR_H_
