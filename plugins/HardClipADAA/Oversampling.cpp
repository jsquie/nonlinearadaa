// Copyright 2024 James Squires
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#include <TargetConditionals.h>
#endif

#include <algorithm>

#include "Oversampling.hpp"
#include "OversamplingStage.hpp"

namespace Oversampling {

Oversampling::~Oversampling() {
  delete[] osBuffer;
  // delete[] up_sample_stages;
  // delete[] down_sample_stages;
  for (auto i = 0; i < factor / 2; ++i) {
    delete up_sample_stages[i];
    delete down_sample_stages[i];
  }
};

void Oversampling::init(const int &newFactor, const int &numSamples,
                        const int &filter_kernel_size,
                        const double *&filter_kernel) {
  factor = newFactor;
  nSamples = numSamples;
  M = filter_kernel_size;
  fKernelBuf = filter_kernel;
  osBuffer = new double[nSamples * factor];
  for (auto i = 0; i < factor / 2; ++i) {
    double *allocated_data_up = new double[M];
    double *allocated_data_down = new double[M];
    up_sample_stages[i] = new OversamplingStage(M, allocated_data_up);
    down_sample_stages[i] = new OversamplingStage(M, allocated_data_down);
  }
};

void Oversampling::processSamplesUp(const double *&input) {
  for (int n = 0; n < nSamples; ++n) {
    for (int j = 0; j < factor / 2; ++j) {
      up_sample_stages[j]->y0 = 0.0f;
      up_sample_stages[j]->y1 = 0.0f;

      up_sample_stages[j]->data[up_sample_stages[j]->pos] =
          static_cast<double>(input[n]);

#if TARGET_OS_MAC

      vDSP_dotprD(fKernelBuf + up_sample_stages[j]->pos, 2,
                  up_sample_stages[j]->data, 1, &up_sample_stages[j]->y0,
                  (up_sample_stages[j]->size - up_sample_stages[j]->pos) / 2);
      vDSP_dotprD(
          fKernelBuf, 2,
          &up_sample_stages[j]
               ->data[up_sample_stages[j]->size - up_sample_stages[j]->pos],
          1, &up_sample_stages[j]->y1, (up_sample_stages[j]->pos) / 2);

#endif

      up_sample_stages[j]->pos = (up_sample_stages[j]->pos == 0)
                                     ? up_sample_stages[j]->size - 1
                                     : up_sample_stages[j]->pos - 1;
      osBuffer[(n << 1) + j] =
          up_sample_stages[j]->y0 + up_sample_stages[j]->y1;
    }
  }
};

void Oversampling::processSamplesDown(float *&output) {
  for (int n = 0; n < nSamples; ++n) {
    // initialize output sample n
    output[n] = 0.0f;

    for (int j = 0; j < factor / 2; ++j) {
      // initialize calcs for CircularBuffer
      down_sample_stages[j]->y0 = 0.0f;
      down_sample_stages[j]->y1 = 0.0f;

      // take every jth sample
      down_sample_stages[j]->data[down_sample_stages[j]->pos] =
          osBuffer[(n << 1) + j];

#if TARGET_OS_MAC

      // dot product every jth sample with every jth kernel value
      // TODO: how does increment change with factor > 2?
      vDSP_dotprD(fKernelBuf + down_sample_stages[j]->pos, 2,
                  down_sample_stages[j]->data, 1, &down_sample_stages[j]->y0,
                  (M - down_sample_stages[j]->pos) / 2);
      vDSP_dotprD(fKernelBuf, 2,
                  &down_sample_stages[j]->data[M - down_sample_stages[j]->pos],
                  1, &down_sample_stages[j]->y1,
                  (down_sample_stages[j]->pos) / 2);

#endif

      // iterate position arguments -- keep it bounded between 0 and size
      down_sample_stages[j]->pos = (down_sample_stages[j]->pos == 0)
                                       ? down_sample_stages[j]->size - 1
                                       : down_sample_stages[j]->pos - 1;

      // increment outbuf[n], because we will have factor number of
      // contributions to outbuf[n]
      output[n] += static_cast<float>(down_sample_stages[j]->y0 +
                                      down_sample_stages[j]->y1);
    }
  }
};

double *Oversampling::getProcessedSamples() { return osBuffer; }

}  // namespace Oversampling
