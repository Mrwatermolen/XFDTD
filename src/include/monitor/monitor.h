#ifndef __MONITOR_H__
#define __MONITOR_H__

#include <memory>

#include "shape/shape.h"
namespace xfdtd {

template <typename T>
class Monitor {
 public:
  Monitor() = default;
  virtual ~Monitor() = default;
  virtual void setData(const T& data);
  virtual void setData(T&& data);
  virtual std::shared_ptr<T> getData();

 protected:
  std::shared_ptr<T> _data;
};

template <typename T>
void Monitor<T>::setData(const T& data) {
  _data = std::make_shared<T>(data);
}

template <typename T>
inline void Monitor<T>::setData(T&& data) {
  _data = std::make_shared<T>(std::move(data));
}

template <typename T>
inline std::shared_ptr<T> Monitor<T>::getData() {
  return _data;
}

template <typename T>
class FieldMonitor : public Monitor<T> {
 public:
  enum class FieldComponent { Ex, Ey, Ez, Hx, Hy, Hz };

  FieldMonitor(Shape shape, FieldComponent field_component);
  virtual ~FieldMonitor() = default;

  std::unique_ptr<Shape> getMonitoredArea() {
    return std::make_unique<Shape>(*_area);
  }

  inline FieldComponent getComponent() { return _component; }

 private:
  std::unique_ptr<Shape> _area;  // monitored areas
  FieldComponent _component;
};

template <typename T>
FieldMonitor<T>::FieldMonitor(Shape shape, FieldComponent field_component)
    : _area(std::make_unique<Shape>(std::move(shape))),
      _component(field_component) {}

}  // namespace xfdtd

#endif
