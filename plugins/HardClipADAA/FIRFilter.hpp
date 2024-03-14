// Copyright 2024 James Squires
//

#pragma once

#include <numeric>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#include <TargetConditionals.h>
#else
#define PAR std::execution::par
#endif

#include <iostream>

#include "./Utils.hpp"

namespace JSCDSP::FIRFilter {

template <typename NumType>
struct CircularBuffer {
  NumType *data;
  unsigned int pos;
  unsigned int size;
};

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

double process_single_sample(CircularBuffer<double> cBuf, double *kernel,
                             unsigned int kernel_size) {
  double result1 = 0.0;
  double result2 = 0.0;

  // std::cout << "Convolving single sample!" << std::endl;
  return result1 + result2;
}

void linear_convolve(const float *input, double *kernel,
                     CircularBuffer<double> circ_buf, double *output,
                     const unsigned int &input_size,
                     const unsigned int &kernel_size) {
  assert(input != nullptr);
  assert(kernel != nullptr);
  assert(output != nullptr);

  // std::cout << "Preparing to convolve signal" << std::endl;
  // std::cout << "circ_buf:" << std::endl;
  // for (auto i = 0u; i < kernel_size; ++i) {
  // std::cout << circ_buf[i] << std::endl;
  // }
  double y = 0.0f;
  double y2 = 0.0f;

  for (auto i = 0u; i < input_size; ++i) {
    circ_buf.data[circ_buf.pos] = input[i];

#if TARGET_OS_MAC

    vDSP_dotprD(kernel + circ_buf.pos, 1, circ_buf.data, 1, &y,
                kernel_size - circ_buf.pos);
    vDSP_dotprD(kernel, 1, &circ_buf.data[kernel_size - circ_buf.pos], 1, &y2,
                circ_buf.pos);

#else

// TODO: THIS IS NOT CORRECT -- needs to be fixed such that it is like above
    std::vector<float> circ_part_1;

    for (int i = static_cast<signed int>(cBuf.pos); i >= 0; --i) {
      circ_part_1.push_back(cBuf.data[i]);
    }

    for (auto i = kernel_size - 1; i > cBuf.pos; --i) {
      circ_part_1.push_back(cBuf.data[i]);
    }

    result1 = std::transform_reduce(PAR circ_part_1.begin(), circ_part_1.end(),
                                    kernel, 0.0f);

#endif

    circ_buf.pos = circ_buf.pos == 0 ? kernel_size - 1 : circ_buf.pos - 1;

    output[i] = y + y2;
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
void buildSinc(NumType *output, int length, NumType cutoff,
               NumType fs) {
  NumType n = 0.0 - static_cast<NumType>(length / 2);

  for (auto i = 0; i < length; ++i) {
    NumType pi_n = Constants::PI * n;
    // std::cout << "pin: " << pi_n << std::endl;
    output[i] = (i == length / 2) ? 1.0 : std::sin( pi_n * cutoff * (fs / 2) ) / pi_n;
    n += 1.0;
  }
}

/**
 * designFIRKaiserKernel takes desired characteristics of
 * FIR window and returns the calculated kernel coefficients
 * in the output buffer
 * \param: output -- size M -- populated with filter coefficients
 * \param: window -- size M -- used for storing window values
 * \param: fs -- current sample rate
 * \param M -- size of kernel
 * \param shape -- beta value of kaiser window
 **/
// results in normalized filter kernel stored in output
template <typename NumType>
void designFIRKaiserKernel(NumType *output, NumType *window, NumType fs,
                           unsigned int M, NumType shape) {
  auto normalized_cutoff = fs / 4;
  // TODO: fix normalized cutoff value
  buildSinc<NumType>(output, M, 0.01, fs);
  // std::cout << "After build sinc" << std::endl;
  // for (auto n = 0u; n < M; ++n) {
    // std::cout << output[n] << std::endl;
  // }
  buildWindow<NumType>(window, M, shape);
  // std::cout << "After build window" << std::endl;
  // for (auto n = 0u; n < M; ++n) {
    // std::cout << window[n] << std::endl;
  // }
  applyWindow<NumType>(output, window, M);

  NumType sum = 0.0;

  // normalize for 1 at DC
  for (auto n = 0u; n < M; ++n) {
    // std::cout << "output before norm: " << output[n] << std::endl;
    sum += output[n];
  }
  for (auto n = 0u; n < M; ++n) {
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
unsigned int computeKaiserLength(NumType width, NumType alpha) {
  return static_cast<unsigned int>(
      (1.0 + (2 * std::sqrt((Constants::PI_SQRD + (alpha * alpha)) /
                            (Constants::PI * width)))));
}

}  // namespace JSCDSP::FIRFilter
