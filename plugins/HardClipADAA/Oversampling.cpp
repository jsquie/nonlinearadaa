// Copyright 2024 James Squires
#include <memory>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#include <TargetConditionals.h>
#endif

#include "FIRFilter.hpp"
#include "Oversampling.hpp"
#include "OversamplingStage.hpp"

constexpr double NEG_ONE_DB = 0.8912509381337456;

namespace Oversampling {

void Oversampling::init(const int &initOSFactor, const int &initNSamples) {
  factor = initOSFactor;
  nSamples = initNSamples;

  up_sample_stage = std::shared_ptr<OversamplingStage>(
      new OversamplingStage(UP_FILTER_TAP_NUM));
  down_sample_stage = std::shared_ptr<OversamplingStage>(
      new OversamplingStage(UP_FILTER_TAP_NUM));

  up_odd_delay_buffer =
      std::shared_ptr<CircularBuffer>(new CircularBuffer(up_odd_delay));
  down_odd_delay_buffer =
      std::shared_ptr<CircularBuffer>(new CircularBuffer(down_odd_delay));

  up_sample_stage->reset();
  down_sample_stage->reset();
}

double Oversampling::delay(const double &input,
                           std::shared_ptr<CircularBuffer> &delay_buf) {
  delay_buf->data[delay_buf->pos] = input * oddidx_scale;
  delay_buf->pos =
      (delay_buf->pos == 0) ? delay_buf->size - 1 : delay_buf->pos - 1;
  return delay_buf->data[delay_buf->pos];
}

double Oversampling::convolve(const double &input,
                              std::shared_ptr<OversamplingStage> &stage,
                              const double kernel[]) {
  stage->y0 = 0.0f;
  stage->y1 = 0.0f;

  const int zPtr = stage->pos;
  const int M = stage->size;

  stage->data[zPtr] = input;

#if TARGET_OS_MAC

  vDSP_dotprD(stage->data.get() + zPtr, 1, kernel, 1, &stage->y0, M - zPtr);
  vDSP_dotprD(stage->data.get(), 1, kernel + (M - zPtr), 1, &stage->y1, zPtr);

#endif

  stage->pos = (stage->pos == 0) ? M - 1 : zPtr - 1;
  return stage->y0 + stage->y1;
}

void Oversampling::processSamplesUp(const float *input,
                                    double *const &osBuffer) {
  assert(input != nullptr);
  // assert(factor == 2);

  for (int n = 0; n < nSamples; ++n) {
    // even indices of input are convolved with up_ftaps and stored in even
    // indices of osBuffer odd indices of input are delayed by
    // up_odd_idx_delay
    for (int j = 0; j < factor; ++j) {
      double y = (j == 0) ? convolve(input[n], up_sample_stage, up_ftaps)
                          : delay(input[n], up_odd_delay_buffer);
      osBuffer[(n << 1) + j] = y * 2;
    }
  }
};

void Oversampling::processSamplesDown(float *const &output,
                                      double *const &osBuffer) {
  assert(output != nullptr);
  assert(osBuffer != nullptr);

  for (int n = 0; n < nSamples; ++n) {
    double y = 0.0;
    for (int j = 0; j < factor; ++j) {
      double inp = osBuffer[(n << 1) + j];
      y += (j == 0) ? convolve(inp, down_sample_stage, up_ftaps)
                    : delay(inp, down_odd_delay_buffer);
    }
    output[n] = y * NEG_ONE_DB;
  }
}
}  // namespace Oversampling
