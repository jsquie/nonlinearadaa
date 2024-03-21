#pragma once
#include <algorithm>
#include <memory>

struct OversamplingStage {
  double y0{0.0f};
  double y1{0.0f};
  std::shared_ptr<double[]> data = {};
  int pos{0};
  const int size{0};

  explicit OversamplingStage(const int &newSize)
      : size(newSize), y0(0.0), y1(0.0), pos(0.0) {
    data = std::shared_ptr<double[]>(new double[newSize]);
  };

  void reset() { std::fill(data.get(), data.get() + size, 0.0f); }
};

struct CircularBuffer {
  std::shared_ptr<double[]> data = {};
  int pos{0};
  const int size{0};

  explicit CircularBuffer(const int &initSize)
      : data(std::shared_ptr<double[]>(new double[initSize])),
        pos(0),
        size(initSize) {
    std::fill(data.get(), data.get() + size, 0.0f);
  };

};
