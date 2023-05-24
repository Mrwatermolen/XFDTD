#ifndef __MONITOR_H__
#define __MONITOR_H__

#include <filesystem>
#include <memory>

#include "electromagnetic_field/electromagnetic_field.h"
#include "shape/shape.h"
#include "util/type_define.h"
namespace xfdtd {

class Monitor {
 public:
  Monitor() = default;
  Monitor(std::unique_ptr<Shape> shape, std::filesystem::path output_dir_path,
          std::string output_file_name);
  Monitor(const Monitor& other) = delete;
  Monitor(Monitor&& other) noexcept = default;
  Monitor& operator=(const Monitor& other) = delete;
  Monitor& operator=(Monitor&& other) noexcept = default;
  virtual ~Monitor() = default;

  virtual const std::unique_ptr<Shape>& getShape() const;
  virtual const std::filesystem::path& getOutputPath() const;
  virtual const std::string& getOutputFileName() const;

  virtual void setOutputDirPath(const std::string& output_dir_path);
  virtual void setOutputFileName(const std::string& output_file_name);

  virtual void setYeeCells(const YeeCellArray& yee_cells) = 0;
  virtual void setYeeCells(YeeCellArray&& yee_cells) = 0;

  virtual void setEMFInstance(std::shared_ptr<EMF> emf) { _emf = std::move(emf); };

  virtual void update(size_t current_time_step) = 0;
  virtual void outputData() = 0;

 protected:
  inline std::shared_ptr<EMF> getEMFInstance() const { return _emf; }

 private:
  std::unique_ptr<Shape> _shape;
  std::filesystem::path _output_dir_path;
  std::string _output_file_name;
  std::shared_ptr<EMF> _emf;
};

}  // namespace xfdtd

#endif
