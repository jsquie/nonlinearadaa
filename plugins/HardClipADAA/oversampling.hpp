// Limited implementation of oversampling

#pragma once

#include <memory>
#include <vector>

#include "OversamplingStage.hpp"

namespace Oversampling {

class Oversampling {
 public:
  void init(const int& initOSFactor, const int& initNSamples, const int& M);

  double convolve(const double& input,
                  std::shared_ptr<OversamplingStage>& stage,
                  const std::shared_ptr<double[]>& kernel);

  void processSamplesUp(const float* input,
                        const std::vector<std::shared_ptr<double[]>>& kernels,
                        double* const& osBuffer);
  void processSamplesDown(float* const& output,
  const std::vector<std::shared_ptr<double[]>>& kernels,
  double* const& osBuffer);

 private:
  std::vector<std::shared_ptr<OversamplingStage>> up_sample_stages;
  std::vector<std::shared_ptr<OversamplingStage>> down_sample_stages;

  int factor;
  int nSamples;
  int M;
};

}  // namespace Oversampling
