// Limited implementation of oversampling

#pragma once

#include <memory>

#include "OversamplingStage.hpp"

namespace Oversampling {

class Oversampling {
 public:
  void init(const int& initOSFactor, const int& initNSamples);
  
  void processSamplesUp(const float* input, double* const cpyBuf, double* const& osBuffer);
  void processSamplesDown(float* const& output, double* const& osBuffer);

 private:
  std::shared_ptr<OversamplingStage> bottom_stage;
  std::vector<OversamplingStage> up_os_stages;
  std::shared_ptr<OversamplingStage> top_stage;
  std::vector<OversamplingStage> down_os_stages;

  int factor;
  int fScale;
  int nSamples;

};

}  // namespace Oversampling
