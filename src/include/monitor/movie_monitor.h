#ifndef _XFDTD_MOVIE_MONITOR_H_
#define _XFDTD_MOVIE_MONITOR_H_

#include <cstddef>
#include <memory>

#include "monitor/monitor.h"
namespace xfdtd {
class MovieMonitor : public Monitor {
 public:
  MovieMonitor(std::unique_ptr<Monitor> monitor, size_t total_time_steps,
               size_t plot_step = 30);

  MovieMonitor(const MovieMonitor& other);

  MovieMonitor& operator=(const MovieMonitor& other);

  MovieMonitor(MovieMonitor&& other) noexcept = default;

  MovieMonitor& operator=(MovieMonitor&& other) noexcept = default;

  ~MovieMonitor() override = default;

  std::unique_ptr<Monitor> clone() const override;

  void init(const std::shared_ptr<const FDTDBasicCoff>& fdtd_basic_coff,
            const std::shared_ptr<const GridSpace>& grid_space,
            const std::shared_ptr<const EMF>& emf) override;

  void update() override;

  void outputData() override;

  const std::string& getOutputPath() const override;

  const std::string& getOutputFileName() const override;

  void setOutputDirPath(const std::string& output_dir_path) override;

  void setOutputFileName(const std::string& output_file_name) override;

  void setPlotStep(size_t plot_step);

 private:
  std::unique_ptr<Monitor> _monitor;
  size_t _total_time_steps;
  size_t _plot_step;
};
}  // namespace xfdtd

#endif  // _XFDTD_MOVIE_MONITOR_H_
