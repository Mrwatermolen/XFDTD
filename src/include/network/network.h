#ifndef _XFDTD_NETWORK_H_
#define _XFDTD_NETWORK_H_

#include <complex>
#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>

#include "network/port.h"

namespace xfdtd {

// only support calculate the calculate port parameters for single excitation
// port
class Network {
 public:
  Network(std::vector<std::unique_ptr<Port>> ports,
          xt::xarray<double> _frequencies, std::string output_path);

  void init(const std::shared_ptr<const FDTDBasicCoff> &fdtd_basic_coff,
            const std::shared_ptr<const GridSpace> &grid_space,
            const std::shared_ptr<const EMF> &emf);

  void outputData();

 private:
  std::vector<std::unique_ptr<Port>> _ports;
  xt::xarray<double> _frequencies;
  std::string _output_path;

  std::unordered_map<size_t, int> _port_map;
  std::unordered_map<int, xt::xarray<std::complex<double>>> _s_parameters;
};
}  // namespace xfdtd

#endif  // _XFDTD_NETWORK_H_
