// Limited implementation of oversampling

#pragma once

#include <memory>

#include "OversamplingStage.hpp"

namespace Oversampling {

class Oversampling {
 public:
  void init(const int& initOSFactor, const int& initNSamples);

  inline double delay(const double& input,
                      std::shared_ptr<CircularBuffer>& delay_buf);

  inline double convolve(const double& input,
                         std::shared_ptr<OversamplingStage>& stage,
                         const double kernel[]);

  void processSamplesUp(const float* input, double* const& osBuffer);
  void processSamplesDown(float* const& output, double* const& osBuffer);

 private:
  std::shared_ptr<OversamplingStage> up_sample_stage;
  std::shared_ptr<OversamplingStage> down_sample_stage;
  std::shared_ptr<CircularBuffer> up_odd_delay_buffer;
  std::shared_ptr<CircularBuffer> down_odd_delay_buffer;

  int factor;
  int nSamples;
};

}  // namespace Oversampling
