#ifndef _MOVIE_MONITOR_H_
#define _MOVIE_MONITOR_H_

#include <cstddef>
#include <filesystem>
#include <memory>

#include "monitor/monitor.h"
namespace xfdtd {
class MovieMonitor : public Monitor {
 public:
  MovieMonitor(std::unique_ptr<Monitor> monitor, size_t total_time_steps,
               size_t plot_step = 30);
  MovieMonitor(const MovieMonitor& other) = delete;
  MovieMonitor(MovieMonitor&& other) noexcept = default;
  MovieMonitor& operator=(const MovieMonitor& other) = delete;
  MovieMonitor& operator=(MovieMonitor&& other) noexcept = default;
  ~MovieMonitor() override = default;

  const std::unique_ptr<Shape>& getShape() const override;
  const std::filesystem::path& getOutputPath() const override;
  const std::string& getOutputFileName() const override;

  void setOutputDirPath(const std::string& output_dir_path) override;
  void setOutputFileName(const std::string& output_file_name) override;
  void setYeeCells(const YeeCellArray& yee_cells) override;
  void setYeeCells(YeeCellArray&& yee_cells) override;

  void update(size_t current_time_step) override;
  void outputData() override;

  inline void setPlotStep(size_t plot_step) { _plot_step = plot_step; }

 private:
  std::unique_ptr<Monitor> _monitor;
  size_t _total_time_steps;
  size_t _plot_step;
};
}  // namespace xfdtd

#endif  // _MOVIE_MONITOR_H_