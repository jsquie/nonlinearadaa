#pragma once 
#include <memory>
#include <algorithm>

struct OversamplingStage {
  double y0{0.0f};
  double y1{0.0f};
  std::shared_ptr<double[]> data;
  int pos{0};
  const int size;

  explicit OversamplingStage(const int &newSize)
      : size(newSize), y0(0.0), y1(0.0), pos(0.0) {
    data = std::shared_ptr<double[]>(new double[newSize]);
  };

  void reset() {
    std::fill(data.get(), data.get() + size, 0.0f);
  }

};

