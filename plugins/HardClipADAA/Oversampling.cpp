// Copyright 2024 James Squires
#include <memory>

#include "SC_InlineUnaryOp.h"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#include <TargetConditionals.h>
#endif

// #include "FIRFilter.hpp"
#include "Oversampling.hpp"
#include "OversamplingStage.hpp"

namespace Oversampling {

void Oversampling::init(const int &initOSFactor, const int &initNSamples,
                        const int &M) {
  factor = initOSFactor;
  nSamples = initNSamples;

  for (int i = 0; i < factor - 1; ++i) {
    up_sample_stages.emplace_back(new OversamplingStage(
        static_cast<int>(std::ceil(static_cast<float>(M) / 2))));
    down_sample_stages.emplace_back(new OversamplingStage(
        static_cast<int>(std::ceil(static_cast<float>(M) / 2))));
  }
  up_sample_stages.emplace_back(
      new OversamplingStage(std::floor(static_cast<double>(M) / 2)));
  down_sample_stages.emplace_back(
      new OversamplingStage(std::floor(static_cast<double>(M) / 2)));

  for (auto st : up_sample_stages) {
    st.reset();
  }

  for (auto st : down_sample_stages) {
    st.reset();
  }
}

double Oversampling::convolve(const double &input,
                              std::shared_ptr<OversamplingStage> &stage,
                              const std::shared_ptr<double[]> &kernel) {
  stage->y0 = 0.0f;
  stage->y1 = 0.0f;

  const int zPtr = stage->pos;
  const int M = stage->size;

  stage->data[zPtr] = input;

#if TARGET_OS_MAC

  vDSP_dotprD(stage->data.get() + zPtr, 1, kernel.get(), 1, &stage->y0,
              M - zPtr);
  vDSP_dotprD(stage->data.get(), 1, kernel.get() + (M - zPtr), 1, &stage->y1,
              zPtr);

#endif

  stage->pos = (stage->pos == 0) ? M - 1 : zPtr - 1;
  return stage->y0 + stage->y1;
}

void Oversampling::processSamplesUp(
    const float *input, const std::vector<std::shared_ptr<double[]>> &kernels,
    double *const &osBuffer) {
  assert(input != nullptr);
  assert(!kernels.empty());
  // assert(factor == 2);

  for (int n = 0; n < nSamples; ++n) {
    for (int j = 0; j < factor; ++j) {
      osBuffer[(n << factor) + j] =
          convolve(static_cast<double>(zapgremlins(input[n])),
                   up_sample_stages.at(j), kernels.at(j)); 
          // FIRFilter::filter_gain();
    }
  }
};

void Oversampling::processSamplesDown(
    float *const &output, const std::vector<std::shared_ptr<double[]>> &kernels,
    double *const &osBuffer) {
  assert(output != nullptr);
  assert(osBuffer != nullptr);
  assert(!kernels.empty());

  for (int n = 0; n < nSamples; ++n) {
    output[n] = 0.0;
    for (int j = 0; j < factor; ++j) {
      output[n] += convolve(osBuffer[(n << 1) + j], down_sample_stages.at(j),
                            kernels.at(j));
    }
    // output[n] *= FIRFilter::filter_gain();
  }
}
}  // namespace Oversampling
