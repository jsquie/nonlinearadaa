#pragma once
#include <algorithm>
#include <memory>

#include "FIRFilter.hpp"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#include <TargetConditionals.h>
#endif

constexpr double NEG_ONE_DB = 0.8912509381337456;

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

  inline void push(const double &val) { data[pos] = val; }
  inline void decrementPos() { pos = (pos == 0) ? size - 1 : pos - 1; }
  inline double delay(const double &val) {
    data[pos] = val;
    pos = (pos == 0) ? size - 1 : pos - 1;
    return data[pos];
  };
};

struct OversamplingStage {
  const int size{0};
  std::shared_ptr<double[]> data = {};
  std::shared_ptr<CircularBuffer> convBuf = {};
  std::shared_ptr<CircularBuffer> delayBuf = {};

  explicit OversamplingStage(const int &targetSize, const int &delay)
      : size(targetSize),
        data(std::shared_ptr<double[]>(new double[targetSize])),
        convBuf(std::shared_ptr<CircularBuffer>(
            new CircularBuffer(FILTER_TAP_NUM))),
        delayBuf(std::shared_ptr<CircularBuffer>(new CircularBuffer(delay))) {
    std::fill(data.get(), data.get() + size, 0.0f);
  };

  inline double convolve(const double &input) {
    double y0 = 0.0;
    double y1 = 0.0;

    const auto buf = convBuf.get();
    const auto bufData = buf->data.get();

    buf->push(input);

    const int filter_size = buf->size;
    const int filter_pos = buf->pos;

#if TARGET_OS_MAC

    vDSP_dotprD(bufData + filter_pos, 1, up_taps, 1, &y0,
                filter_size - filter_pos);
    vDSP_dotprD(bufData, 1, up_taps + (filter_size - filter_pos), 1,
                &y1, filter_pos);
#else
    y0 = std::inner_product(buf->data.get() + filter_pos,
                            buf->data.get() + filter_size, up_taps, 0.0);
    y0 = std::inner_product(buf->data.get(), buf->data.get() + filter_pos,
                            up_taps + (filter_size - filter_pos), y0);
#endif
    buf->decrementPos();
    return y0 + y1;
  }

  void processUp(const double *input, const int &input_size) {
    for (int n = 0; n < input_size; ++n) {
      data[n << 1] = convolve(input[n]) * upScaleCoef;
      data[(n << 1) + 1] = delayBuf->delay(input[n] * foldScaleCoef * upScaleCoef);
    }
  };

  void processDown(const double *input, const int &input_size) {
    for (int n = 0; n < input_size; n += 2) {
      double res = 0.0;
      res += convolve(input[n]);
      res += delayBuf->delay(input[n + 1] * foldScaleCoef);
      data[n >> 1] = res;
    }
  }
};
