// PluginHardClipADAA.cpp
// James Squires (squires.jr@gmail.com)

#include "SC_PlugIn.hpp"
#include "HardClipADAA.hpp"
#include <cstring>
#include <assert.h>
#include <iostream>

static InterfaceTable* ft;

namespace HardClipADAA {

HardClipADAA::HardClipADAA() {
  coefs_buffer = (float*)RTAlloc(mWorld, (nTaps + 1) * sizeof(float));
  s1_prevs_buf = (float*)RTAlloc(mWorld, (nTaps + 1) * sizeof(float));
  s2_prevs_buf = (float*)RTAlloc(mWorld, (nTaps + 1) * sizeof(float));
  up_sample_buffer = (float*)RTAlloc(mWorld, inBufferSize(0) * 2 * sizeof(float));
  us_stage1filtered_buf = (float*)RTAlloc(mWorld, inBufferSize(0) * 2 * sizeof(float));
  us_stage2filtered_buf = (float*)RTAlloc(mWorld, inBufferSize(0) * 2 * sizeof(float));
  us_process_buf = (float*)RTAlloc(mWorld, inBufferSize(0) * 2 * sizeof(float));
  set_calc_function<HardClipADAA, &HardClipADAA::next>();


  assert(coefs_buffer != NULL);
  assert(s1_prevs_buf != NULL);
  assert(s2_prevs_buf != NULL);
  assert(up_sample_buffer != NULL);
  assert(us_stage2filtered_buf != NULL);
  assert(us_stage1filtered_buf != NULL);
  assert(us_process_buf != NULL);

  const int fc = 23500;

  // determine coefficient values
  coefs(coefs_buffer, fc);

  s1_prevs_ptr = 0;
  s2_prevs_ptr = 0;

  // zero the prevs buf to start
  memset(s1_prevs_buf, 0.0, (nTaps + 1) * sizeof(float));
  memset(s2_prevs_buf, 0.0, (nTaps + 1) * sizeof(float));

  std::cout << "Instantiating hard clip" << std::endl;
}

HardClipADAA::~HardClipADAA() {
  RTFree(mWorld, coefs_buffer);
  RTFree(mWorld, s1_prevs_buf);
  RTFree(mWorld, s2_prevs_buf);
  RTFree(mWorld, up_sample_buffer);
  RTFree(mWorld, us_stage1filtered_buf);
  RTFree(mWorld, us_stage2filtered_buf);
  RTFree(mWorld, us_process_buf);
}

inline void HardClipADAA::coefs(float* c_buff, const int fc) {
  const float normed_cutoff = fc / (sampleRate() * 2.0);

  for (int k = 0; k <= nTaps; ++k) {
    float k_coef = k - (nTaps * 0.5);
    float twopik = twopi * k_coef;
    float pik = pi * k_coef;
    float lhs = std::sin(twopik * normed_cutoff) / pik; 
    float rhs = 0.54 + (0.46 * std::cos(twopik / nTaps));

    if (k == nTaps / 2) {
      c_buff[k] = 2 * normed_cutoff;
    } else {
      c_buff[k] = lhs * rhs;
    }
  }
}

inline void HardClipADAA::stage_1_filter(int nSamples) {
  const float* input = in(0);
  // filter at 23500 hz for sampleRate * 2 
  for (int i = 0; i < nSamples * 2; ++i) {
    // a_0 * x(n)

    float sample;
    // if idx is even, use input
    // otherwise add 0.0
    if (i % 2 == 0) {
      sample = input[i >> 1];
    } else {
      sample = 0.0;
    }

    float a0_samp = sample * coefs_buffer[0];

    float coef_cals = 0.0;
    // calculate coefs[1..] * prevs[]
    for (int j = 0; j < nTaps; ++j) {
      coef_cals += (coefs_buffer[j + 1] * (s1_prevs_buf[(s1_prevs_ptr + j) % (nTaps + 1)]));
    }

    us_stage1filtered_buf[i] = coef_cals + a0_samp;

    if (s1_prevs_ptr == 0) {
      s1_prevs_ptr = nTaps;
    } else {
      s1_prevs_ptr--;
    }

    s1_prevs_buf[s1_prevs_ptr] = sample;
  }
}

inline void HardClipADAA::stage_2_filter(int nSamples) {
  // filter at 23500 hz for sampleRate * 2 
  for (int i = 0; i < nSamples * 2; ++i) {
    // a_0 * x(n)
    float sample = us_process_buf[i];
    float a0_samp = sample * coefs_buffer[0];

    float coef_cals = 0.0;
    // calculate coefs[1..] * prevs[]
    for (int j = 0; j < nTaps; ++j) {
      coef_cals += (coefs_buffer[j + 1] * (s2_prevs_buf[(s2_prevs_ptr + j) % (nTaps + 1)]));
    }

    us_stage2filtered_buf[i] = coef_cals + a0_samp;

    if (s2_prevs_ptr == 0) {
      s2_prevs_ptr = nTaps;
    } else {
      s2_prevs_ptr--;
    }

    s2_prevs_buf[s2_prevs_ptr] = sample;
  }
}


void HardClipADAA::next(int nSamples) {

  const float* input = in(0);
  // Output buffer
  float* outbuf = out(0);

  // stage_1_filter(nSamples);

  // non linear process
  for (int i = 0; i < nSamples; ++i) {
    float sample = static_cast<double>(input[i]);
    double res = std::tanh(sample);

    outbuf[i] = static_cast<float>(res);
  }
  // stage 2 filter
  // stage_2_filter(nSamples);

  // down sample the result of the processing
  // for (int i = 0; i < nSamples; ++i) {
    // outbuf[i] = us_stage2filtered_buf[i << 1];
  // }

}

} // namespace HardClipADAA

PluginLoad(HardClipADAAUGens) {
  // Plugin magic
  ft = inTable;
  registerUnit<HardClipADAA::HardClipADAA>(ft, "HardClipADAA", true);
}
