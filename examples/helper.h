#ifndef __XFDTD_EXAMPLE_HELPER_H__
#define __XFDTD_EXAMPLE_HELPER_H__

#include <chrono>

namespace xfdtd_example {

template <typename Func>
auto timeSomething(Func &&f) {
  auto start = std::chrono::high_resolution_clock::now();
  f();
  auto end = std::chrono::high_resolution_clock::now();
  return end - start;
}

}  // namespace xfdtd_example

#endif  // __XFDTD_EXAMPLE_HELPER_H__
