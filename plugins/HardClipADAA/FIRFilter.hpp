// Copyright 2024 James Squires
#pragma once

#include "./Utils.hpp"

namespace FIRFilter {

class FIRFilter {

public:
  FIRFilter() = {}



template <typename NumType>
inline NumType zeroethOrderBessel(NumType x) {
  const NumType eps = 0.00001;

  NumType besselValue = 0;
  NumType term = 1;
  NumType m = 0;

  while (term > eps * besselValue) {
    besselValue += term;
    ++m;
    term *= (x * x) / (4 * m * m);
  }

  return besselValue;
}

// win is output buffer to write to, must be of size: sizeof(NumType) * M
// M is window size
// shape is alpha beta parameter determined from computeShape()
template <typename NumType>
void buildWindow(NumType *win, unsigned int M, NumType shape) {
  const NumType oneOverDenom = 1.0 / zeroethOrderBessel<NumType>(shape);
  const unsigned int N = M - 1;
  const NumType oneOverN = 1.0 / N;

  for (unsigned int n = 0; n <= N; ++n) {
    const NumType K = (2.0 * n * oneOverN) - 1.0;
    const double arg = std::sqrt(1.0 - (K * K));

    win[n] = zeroethOrderBessel<NumType>(shape * arg) * oneOverDenom;
  }
}

template <typename NumType>
void applyWindow(NumType *output, NumType *window, unsigned int length) {
  for (unsigned int n = 0; n < length; n++) {
    output[n] = output[n] * window[n];
  }
}

template <typename NumType>
void buildSinc(NumType *output, unsigned int length, NumType cutoff,
               NumType fs) {
  for (int n = 0; n < length; ++n) {
    const NumType scaled_n = n - ((length - 1) / 2);
    output[n] = std::sin(Constants::TWO_PI * cutoff * (scaled_n / fs)) /
                (Constants::TWO_PI * scaled_n);
  }
}

// results in normalized filter kernel stored in output
template <typename NumType>
void designFIRKaiserKernel(NumType *output, NumType *window, NumType fs,
                           unsigned int M, NumType shape) {
  auto normalized_cutoff = fs / 4;
  buildSinc<NumType>(output, M, normalized_cutoff, fs);
  buildWindow<NumType>(window, M, shape);
  applyWindow<NumType>(output, window, M);

  NumType sum = 0.0;

  // normalize for 1 at DC
  for (unsigned int n = 0; n < M; ++n) {
    sum += output[n];
  }
  for (unsigned int n = 0; n < M; ++n) {
    output[n] /= sum;
  }
}

// compute beta given desired passband attenuation
template <typename NumType>
NumType computeKaiserBeta(NumType passband_attenuation) {
  NumType alpha;

  if (passband_attenuation > 60.0) {
    alpha = 0.12438 * (passband_attenuation + 6.3);
  } else if (passband_attenuation > 13.26) {
    alpha = 0.76609 * (std::pow((passband_attenuation - 13.26), 0.4)) +
            0.09834 * (passband_attenuation - 13.26);
  } else {
    alpha = 0.0;
  }

  return alpha;
}

// Compute M (length of window)
// given normalized width (between 0 and 0.5)
// and alpha (desired passband attenuation)
template <typename NumType>
NumType computeKaiserLength(NumType width, NumType alpha) {
  return (1.0 + (2 * std::sqrt((Constants::PI_SQRD + (alpha * alpha)) /
                               (Constants::PI * width))));
}

// convolves input of buffer_size with filter kernel of kernel_size
// places results in output
template <typename NumType>
void filterConvolve(NumType *input, NumType *kernel, NumType *output,
                    unsigned int buffer_size, unsigned int kernel_size, int numSamples) {
  for (int i = 0; i < numSamples, ++i) {
    convolveSingleSample(input[i], , NumType *kernel, unsigned int kernel_size, unsigned int buffer_idx)
  }
}

template <typename NumType>
NumType convolveSingleSample(NumType sample, NumType *circular_buffer, NumType *kernel, unsigned int kernel_size, unsigned int buffer_idx) {

}
}

}  // namespace FIRFilter
