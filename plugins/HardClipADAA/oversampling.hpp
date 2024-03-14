// Limited implementation of oversampling

#pragma once

#include "OversamplingStage.hpp"


namespace Oversampling {


class Oversampling {
 public:

  ~Oversampling();

  void init(const int& newFactor, const int& numSamples, const int& filter_kernel_size, const double*&filter_kernel);

  void processSamplesUp(const double* &input);

  void processSamplesDown(float* &output);

  double* getProcessedSamples();


 private:

  OversamplingStage* up_sample_stages[8];
  OversamplingStage* down_sample_stages[8];

  int factor;
  int nSamples;
  int M;

  double* osBuffer;
  const double* fKernelBuf;

};


} // namespace JSCDSP::Oversampling
