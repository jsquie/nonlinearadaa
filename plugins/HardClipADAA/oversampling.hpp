// Limited implementation of oversampling

#pragma once

#include <vector>

#include "OversamplingStage.hpp"

namespace Oversampling {

class Oversampling {
 public:
  void init(const int& initOSFactor, const int& initNSamples, const int& M);

  double convolve(const double& input, OversamplingStage& stage,
                  std::shared_ptr<double[]> const& kernel);

  void processSamplesUp(const float* const& input,
                        std::vector<std::shared_ptr<double[]>>& kernels,
                        std::vector<double>&osBuffer);

  void processSamplesDown(float* const& output,
                          std::vector<std::shared_ptr<double[]>>& kernels,
                          double* const& osBuffer);

 private:
  std::vector<OversamplingStage> up_sample_stages;
  std::vector<OversamplingStage> down_sample_stages;

  int factor;
  int nSamples;
  int M;
};

}  // namespace Oversampling
