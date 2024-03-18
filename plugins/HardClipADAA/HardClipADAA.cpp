// Copyright 2024 James Squires
// Nonlinear shapers with Anti Derivative Anti Aliasing
// James Squires (squires.jr@gmail.com)
//
// Adapted from Chow DSP
#include "HardClipADAA.hpp"

#include <cmath>
#include <iostream>
#include <memory>

#include "Li2.hpp"
#include "SC_InterfaceTable.h"
#include "SC_PlugIn.hpp"
#include "Utils.hpp"

static InterfaceTable *ft;

static unsigned int numCalls = 0;

namespace JSCDSP::ADAA {

constexpr double TOL = 1.0e-5;

template <typename Func, typename Ad_Func>
static inline double next_first_adaa(const double &s, double &x1,
                                     double &ad1_x1, Func f, Ad_Func first_ad) {
  double res;
  double ad1_x = first_ad(s);
  double diff = s - x1;

  if (std::abs(diff) <= TOL) {
    res = f(0.5 * (s + x1));
  } else {
    res = (ad1_x - ad1_x1) / diff;
  }

  ad1_x1 = ad1_x;
  x1 = s;

  return res;
}

template <typename FuncFirstAD, typename FuncSecondAD>
inline double calcD(const double &v0, const double &x1, double &ad2_x0,
                    const double &ad2_x1, FuncFirstAD f_first_ad,
                    FuncSecondAD f_second_ad) {
  ad2_x0 = f_second_ad(v0);

  if (std::abs(v0 - x1) <= TOL) {
    return f_first_ad(0.5 * (v0 + x1));
  } else {
    return (ad2_x0 - ad2_x1) / (v0 - x1);
  }
}

template <typename Func, typename FuncFirstAD, typename FuncSecondAD>
inline double fallback(const double &x, const double &x1, const double &x2,
                       const double &ad2_x1, Func f, FuncFirstAD f_first_ad,
                       FuncSecondAD f_second_ad) {
  const double xbar = 0.5 * (x + x2);
  const double delta = xbar - x1;

  if (std::abs(delta) <= TOL) {
    return f(0.5 * (xbar + x1));
  } else {
    return (2.0 / delta) *
           (f_first_ad(xbar) + ((ad2_x1 - f_second_ad(xbar)) / delta));
  }
}

template <typename Func, typename FuncFirstAD, typename FuncSecondAD>
inline double next_second_adaa(const double &s, double &x1, double &x2,
                               double &ad2_x0, double &ad2_x1, double &d2,
                               Func f, FuncFirstAD f_first_ad,
                               FuncSecondAD f_second_ad) {
  double res;
  double d1 = calcD(s, x1, ad2_x0, ad2_x1, f_first_ad, f_second_ad);

  // x_n - x_{n-1} <= epsilon, then use
  if (std::abs(s - x2) <= TOL) {  // fallback
    res = fallback(s, x1, x2, ad2_x1, f, f_first_ad, f_second_ad);
  } else {
    res = (2.0 / (s - x2)) * (d1 - d2);
  }

  d2 = d1;
  x2 = x1;
  x1 = s;
  ad2_x1 = ad2_x0;

  return res;
}

enum AntiDerivativeLevel { First = 1, Second = 2 };

}  // namespace JSCDSP::ADAA

namespace JSCDSP::HardClipADAA {

HardClipADAA::HardClipADAA() {
  M = 63;
  osBuffer = reinterpret_cast<double *>(
      RTAlloc(mWorld, inBufferSize(0) * in0(OverSample) * sizeof(double)));

  float Mf = static_cast<float>(M);
  int first_buffer_size = std::ceil(Mf / in0(OverSample));
  int last_buffer_size = std::floor(Mf / in0(OverSample));

/*

FIR filter designed with
http://t-filter.appspot.com

sampling frequency: 44100 Hz

* 0 Hz - 1000 Hz
  gain = 1
  desired ripple = 5 dB
  actual ripple = 4.050040059095146 dB

* 1500 Hz - 22050 Hz
  gain = 0
  desired attenuation = -25 dB
  actual attenuation = -25.256251343510396 dB

*/

#define FILTER_TAP_NUM 63

static double filter_taps[FILTER_TAP_NUM] = {
  -0.030660379541817155,
  -0.006930252652193655,
  -0.007426606875288303,
  -0.0077304521574798925,
  -0.007805654047780043,
  -0.007627358073715853,
  -0.007123361170733022,
  -0.006377495207966196,
  -0.005212485606290447,
  -0.0037179612791280836,
  -0.0018893069905069938,
  0.0003031347728881851,
  0.0028518138446229848,
  0.005755694749310335,
  0.008945405753532943,
  0.012416984304270736,
  0.016119314979375885,
  0.019997655165236015,
  0.02401203285079631,
  0.02809885339907276,
  0.03220543281661772,
  0.03624236586185808,
  0.040151565039062635,
  0.043869271258687656,
  0.047323200707502136,
  0.050460686719297834,
  0.05321950592837737,
  0.05555898793967278,
  0.05742296873168411,
  0.05877846833665994,
  0.059606627462459476,
  0.059882326502100264,
  0.059606627462459476,
  0.05877846833665994,
  0.05742296873168411,
  0.05555898793967278,
  0.05321950592837737,
  0.050460686719297834,
  0.047323200707502136,
  0.043869271258687656,
  0.040151565039062635,
  0.03624236586185808,
  0.03220543281661772,
  0.02809885339907276,
  0.02401203285079631,
  0.019997655165236015,
  0.016119314979375885,
  0.012416984304270736,
  0.008945405753532943,
  0.005755694749310335,
  0.0028518138446229848,
  0.0003031347728881851,
  -0.0018893069905069938,
  -0.0037179612791280836,
  -0.005212485606290447,
  -0.006377495207966196,
  -0.007123361170733022,
  -0.007627358073715853,
  -0.007805654047780043,
  -0.0077304521574798925,
  -0.007426606875288303,
  -0.006930252652193655,
  -0.030660379541817155
};

  // os.init(in0(OverSample), inBufferSize(0), M);
  for (int i = 0; i < in0(OverSample); ++i) {
    if (i == in0(OverSample) - 1) {
      kernels.emplace_back(new double[last_buffer_size]);
      for (int j = 0; j < last_buffer_size; ++j) {
        kernels.at(i).get()[j] = filter_taps[(j << 1) + i];
      }
    } else {
      kernels.emplace_back(new double[first_buffer_size]);
      for (int j = 0; j < first_buffer_size; ++j) {
        kernels.at(i).get()[j] = filter_taps[(j << 1) + i];
      }
    }
  }

  os.init(in0(OverSample), inBufferSize(0), M);

  x1 = 0.0;
  x2 = 0.0;
  d2 = 0.0;
  ad1_x1 = 0.0;
  ad2_x1 = 0.0;
  ad2_x0 = 0.0;

  mCalcFunc = make_calc_function<HardClipADAA, &HardClipADAA::next_aa>();
  next_aa(1);
}

HardClipADAA::~HardClipADAA() { RTFree(mWorld, osBuffer); }

inline double HardClipADAA::clip(const double &v) {
  return (std::abs(v) <= 1.0f) ? v : 1.0;
}

inline double HardClipADAA::hc_first_ad(const double &v) {
  return (std::abs(v) <= 1.0f) ? (v * v) * 0.5 : (v * signum(v)) - 0.5;
}

inline double HardClipADAA::hc_second_ad(const double &v) {
  if (std::abs(v) <= 1.0f) {
    return std::pow(v, 3) / 6.0;
  } else {
    return (((std::pow(v, 2) * 0.5) + (1.0 / 6.0)) * signum(v)) - (v * 0.5);
  }
}

void HardClipADAA::next_aa(int nSamples) {
  int twoNSamples = nSamples * 2;

  const float *sig = in(Input);
  float *outbuf = out(Out1);
  const int adlevel = static_cast<int>(in0(AntiDerivativeLevel));

  std::vector<double> osB;

  // up sample -- results are stored in osBuffer
  os.processSamplesUp(sig, kernels, osB);

  // assert(osB.size() == 128);
  // process
  for (auto it = osB.begin(); it != osB.end(); ++it) {
    if (adlevel == JSCDSP::ADAA::AntiDerivativeLevel::First) {
      *it = JSCDSP::ADAA::next_first_adaa(*it, x1, ad1_x1,
                                                  HardClipADAA::clip,
                                                  HardClipADAA::hc_first_ad);

    } else if (adlevel == ADAA::AntiDerivativeLevel::Second) {
      *it = JSCDSP::ADAA::next_second_adaa(
          *it , x1, x2, ad2_x0, ad2_x1, d2, HardClipADAA::clip,
          HardClipADAA::hc_first_ad, HardClipADAA::hc_second_ad);
    }
  }

  int i = 0;
  for (auto it = osB.cbegin(); it != osB.cend(); it += 2) {
    outbuf[i] = static_cast<float>(*it); 
    ++i;
  }

  // os.processSamplesDown(outbuf, kernels, osBuffer);
}

}  // namespace JSCDSP::HardClipADAA

namespace JSCDSP::TanhADAA {

TanhADAA::TanhADAA() {
  x1 = 0.0;
  x2 = 0.0;
  d2 = 0.0;
  ad1_x1 = 0.0;
  ad2_x1 = 0.0;
  ad2_x0 = 0.0;

  mCalcFunc = make_calc_function<TanhADAA, &TanhADAA::next_aa>();
  next_aa(1);
}

TanhADAA::~TanhADAA() {}

inline double TanhADAA::tanh_first_ad(const double &v) {
  return std::log(std::cosh(v));
}

inline double TanhADAA::tanh_second_ad(const double &v) {
  const auto exp = std::exp(-2 * v);
  return 0.5 * (static_cast<double>(polylogarithm::Li2(-exp)) -
                v * (v + 2.0 * std::log(exp + 1.) -
                     2.0 * std::log(std::cosh(v)))) +
         (Constants::PI_SQRD / 24.0);
}

void TanhADAA::next_aa(int nSamples) {
  const float *sig = in(Input);

  float *outbuf = (float *)out(Out1);
  const int adlevel = static_cast<int>(in0(AntiDerivativeLevel));

  /**
  os.processSamplesUp();

  // process
  for (auto s : os.getProcessedSamples()) {
    if (adlevel == JSCDSP::ADAA::AntiDerivativeLevel::First) {
      s = JSCDSP::ADAA::next_first_adaa(
          s, &x1, &ad1_x1, [](double v) { return std::tanh(v); },
          TanhADAA::tanh_first_ad);

    } else if (adlevel == ADAA::AntiDerivativeLevel::Second) {
      s = JSCDSP::ADAA::next_second_adaa(
          s, &x1, &x2, &ad2_x0, &ad2_x1, &d2,
          [](double v) { return std::tanh(v); }, TanhADAA::tanh_first_ad,
          TanhADAA::tanh_second_ad);
    }
  }

  os.processSamplesDown();

  auto down_sampled = os.getProcessedSamples();
  // downsample and out
  for (int i = 0; i < nSamples; ++i) {
    outbuf[i] = static_cast<float>(down_sampled[i]);
  }
  **/
}

}  // namespace JSCDSP::TanhADAA

PluginLoad(HardClipADAAUGens) {
  // Plugin magic
  ft = inTable;
  registerUnit<JSCDSP::HardClipADAA::HardClipADAA>(ft, "HardClipADAA", false);
  registerUnit<JSCDSP::TanhADAA::TanhADAA>(ft, "TanhADAA", false);
}
