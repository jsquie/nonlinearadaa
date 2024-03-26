// Copyright 2024 James Squires
#include "Oversampling.hpp"

#include <algorithm>
#include <memory>

#include "FIRFilter.hpp"
#include "OversamplingStage.hpp"


namespace Oversampling {

void Oversampling::init(const int &initOSFactor, const int &initNSamples) {
  factor = initOSFactor;
  fScale = std::pow(2, factor);
  nSamples = initNSamples;
  bottom_stage = std::shared_ptr<OversamplingStage>(
      new OversamplingStage(nSamples * 2, up_delay));
  top_stage = std::shared_ptr<OversamplingStage>(
      new OversamplingStage((nSamples * fScale) / 2, down_delay));

  for (auto i = 1; i < factor; ++i) {
    up_os_stages.emplace_back(nSamples * std::pow(2, i + 1), up_delay);
    down_os_stages.emplace_back(nSamples * std::pow(2, factor - (i + 1)),
                                down_delay);
  }
}

void Oversampling::processSamplesUp(const float *input, double *cpyBuf,
                                    double *const &osBuffer) {
  auto last_stage = bottom_stage.get();
  std::transform(input, input + nSamples, cpyBuf,
                 [](const float &val) { return static_cast<double>(val); });
  last_stage->processUp(cpyBuf, nSamples);
  for (auto &stage : up_os_stages) {
    stage.processUp(last_stage->data.get(), last_stage->size);
    last_stage = &stage;
  }
  std::copy_n(last_stage->data.get(), nSamples * fScale, osBuffer);
};

void Oversampling::processSamplesDown(float *const &output,
                                      double *const &osBuffer) {
  auto last_stage = top_stage.get();
  last_stage->processDown(osBuffer, nSamples * fScale);
  for (auto &stage : down_os_stages) {
    stage.processDown(last_stage->data.get(), last_stage->size);
    last_stage = &stage;
  }
  for (int i = 0; i < nSamples; ++i) {
    output[i] = static_cast<float>(last_stage->data[i]);
  }
};

}  // namespace Oversampling
