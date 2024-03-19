// Copyright 2024 James Squires
// Nonlinear shapers with Anti Derivative Anti Aliasing
// James Squires (squires.jr@gmail.com)
//
// Adapted from Chow DSP
#include "HardClipADAA.hpp"

#include <cmath>

#include "FIRFilter.hpp"
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

  // M = FIRFilter::FILTER_TAP_NUM;
  M = 209;


  osBuffer = reinterpret_cast<double *>(
      RTAlloc(mWorld, inBufferSize(0) * in0(OverSample) * sizeof(double)));

  float Mf = static_cast<float>(M);
  int first_buffer_size = std::ceil(Mf / in0(OverSample));
  int last_buffer_size = std::floor(Mf / in0(OverSample));

  for (int i = 0; i < in0(OverSample); ++i) {
    if (i == in0(OverSample) - 1) {
      kernels.emplace_back(new double[last_buffer_size]);
      for (int j = 0; j < last_buffer_size; ++j) {
        kernels.at(i).get()[j] = FIRFilter::filter_taps[(j << 1) + i];
      }
    } else {
      kernels.emplace_back(new double[first_buffer_size]);
      for (int j = 0; j < first_buffer_size; ++j) {
        kernels.at(i).get()[j] = FIRFilter::filter_taps[(j << 1) + i];
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

  // up sample -- results are stored in osBuffer
  os.processSamplesUp(sig, kernels, osBuffer);

  // assert(osB.size() == 128);
  // process
  for (auto i = 0; i < twoNSamples; ++i) {
    if (adlevel == JSCDSP::ADAA::AntiDerivativeLevel::First) {
      osBuffer[i] = JSCDSP::ADAA::next_first_adaa(osBuffer[i], x1, ad1_x1,
                                                  HardClipADAA::clip,
                                                  HardClipADAA::hc_first_ad);

    } else if (adlevel == ADAA::AntiDerivativeLevel::Second) {
      osBuffer[i] = JSCDSP::ADAA::next_second_adaa(
          osBuffer[i], x1, x2, ad2_x0, ad2_x1, d2, HardClipADAA::clip,
          HardClipADAA::hc_first_ad, HardClipADAA::hc_second_ad);
    }
  }

  os.processSamplesDown(outbuf, kernels, osBuffer);
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
