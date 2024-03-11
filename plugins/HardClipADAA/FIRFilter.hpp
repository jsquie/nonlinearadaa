// Copyright 2024 James Squires
//

#pragma once

#include <exception>
#define PAR std::execution::par

#ifdef __APPLE__
#include <TargetConditionals.h>
#include <Accelerate/Accelerate.h>
#endif

#include "./Utils.hpp"

namespace JSCDSP::FIRFilter {

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

double process_single_sample(double *circular_buf, double *kernel,
                                     unsigned int kernel_size,
                                     unsigned int pos) {
  double result1 = 0.0;
  double result2 = 0.0;

#if TARGET_OS_MAC

  vDSP_dotprD(&circular_buf[pos], -1, kernel, 1, &result1, pos + 1);
  vDSP_dotprD(&circular_buf[kernel_size - 1], -1, &kernel[pos+1], 1, &result2, kernel_size - pos - 1);

#else

  std::vector<float> circ_part_1;

  for (int i = static_cast<signed int>(pos); i >= 0; --i) {
    circ_part_1.push_back(circ_buf[i]);
  }

  for (auto i = kernel_size - 1; i > pos; --i) {
    circ_part_1.push_back(circ_buf[i]);
  }

  result1 = std::transform_reduce(PAR circ_part_1.begin(), circ_part_1.end(), kernel, 0.0f);

#endif

  return result1 + result2;
}

void linear_convolve(float *input, float *kernel, float* circ_buf, float* output,
                                   const unsigned int &input_size,
                                   const unsigned int &kernel_size) {

  unsigned int pos = 0;

  for (unsigned int i = 0; i < input_size + kernel_size - 1; ++i) {

    circ_buf[pos] = (i < input_size) ? input[i] : 0.0;

    output[i] = process_single_sample(circ_buf, kernel, kernel_size, pos);

    pos = (++pos) % kernel_size;
  }

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

}  // namespace JSCDSP::FIRFilter
