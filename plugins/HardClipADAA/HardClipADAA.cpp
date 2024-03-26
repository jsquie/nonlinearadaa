// Copyright 2024 James Squires
// Nonlinear shapers with Anti Derivative Anti Aliasing
// James Squires (squires.jr@gmail.com)
//
// Adapted from Chow DSP
#include "HardClipADAA.hpp"

#include <cmath>

#include "Li2.hpp"
#include "SC_InterfaceTable.h"
#include "SC_PlugIn.hpp"
#include "Utils.hpp"

static InterfaceTable *ft;
#define print(x) std::cout << x << ", "
#define printEln() std::cout << "\n"
#define printfl() std::cout << "" std::endl
static unsigned int numCalls = 0;

namespace JSCDSP::ADAA {

constexpr double TOL = 1.0e-5;

template <typename Func, typename Ad_Func>
inline double next_first_adaa(const double &s, double &x1, double &ad1_x1,
                              Func f, Ad_Func first_ad) {
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
  // float M = static_cast<float>(M);
  const int osFactor = in0(OverSample);
  const int osScale = std::pow(2, osFactor);
  fScale = osScale;
  osBuffer = static_cast<double *>(
      RTAlloc(mWorld, inBufferSize(0) * osScale * sizeof(double)));

  cpyBuf =
      static_cast<double *>(RTAlloc(mWorld, inBufferSize(0) * sizeof(double)));

  os.init(osFactor, inBufferSize(0));

  x1 = 0.0;
  x2 = 0.0;
  d2 = 0.0;
  ad1_x1 = 0.0;
  ad2_x1 = 0.0;
  ad2_x0 = 0.0;

  mCalcFunc = make_calc_function<HardClipADAA, &HardClipADAA::next_aa>();
  next_aa(1);
}

HardClipADAA::~HardClipADAA() {
  RTFree(mWorld, osBuffer);
  RTFree(mWorld, cpyBuf);
}

inline double HardClipADAA::signum(const double &d) {
  return (d < 0.0) ? -1.0 : 1.0;
}

inline double HardClipADAA::clip(const double &d) {
  const double t = d < -1.0 ? -1.0 : d;
  return t > 1.0 ? 1.0 : t;
}

inline double HardClipADAA::hc_first_ad(const double &v) {
  return (std::abs(v) <= 1.0) ? (v * v) * 0.5 : (v * signum(v)) - 0.5;
}

inline double HardClipADAA::hc_second_ad(const double &v) {
  const double oneSixth = 1.0 / 6.0;
  if (std::abs(v) <= 1.0) {
    return (v * v * v) * oneSixth;
  } else {
    return (((v * v * 0.5) + oneSixth) * signum(v)) - (v * 0.5);
  }
}

void HardClipADAA::next_aa(int nSamples) {
  // println("Starting a next_aa call");
  const float *sig = in(Input);
  float *outbuf = out(Out1);
  const int adlevel = static_cast<int>(in0(AntiDerivativeLevel));
  const int osBufferSize = nSamples * fScale;

  // up sample -- results are stored in osBuffer
  // up sample -- results are stored in osBuffer
  os.processSamplesUp(sig, cpyBuf, osBuffer);

  // process
  for (auto i = 0; i < osBufferSize; ++i) {
    if (adlevel == ADAA::AntiDerivativeLevel::First) {
      double curr =
          ADAA::next_first_adaa(osBuffer[i], x1, ad1_x1, HardClipADAA::clip,
                                HardClipADAA::hc_first_ad);
      osBuffer[i] = curr;

    } else if (adlevel == ADAA::AntiDerivativeLevel::Second) {
      double curr = ADAA::next_second_adaa(
          osBuffer[i], x1, x2, ad2_x0, ad2_x1, d2, HardClipADAA::clip,
          HardClipADAA::hc_first_ad, HardClipADAA::hc_second_ad);
      osBuffer[i] = curr;
    }
  }

  os.processSamplesDown(outbuf, osBuffer);
}

}  // namespace JSCDSP::HardClipADAA

namespace JSCDSP::TanhADAA {

TanhADAA::TanhADAA() {
  // float M = static_cast<float>(M);
  const int osFactor = in0(OverSample);
  const int osScale = std::pow(2, osFactor);
  fScale = osScale;
  osBuffer = static_cast<double *>(
      RTAlloc(mWorld, inBufferSize(0) * osScale * sizeof(double)));

  cpyBuf =
      static_cast<double *>(RTAlloc(mWorld, inBufferSize(0) * sizeof(double)));

  os.init(osFactor, inBufferSize(0));

  x1 = 0.0;
  x2 = 0.0;
  d2 = 0.0;
  ad1_x1 = 0.0;
  ad2_x1 = 0.0;
  ad2_x0 = 0.0;

  mCalcFunc = make_calc_function<TanhADAA, &TanhADAA::next_aa>();
  next_aa(1);
}

TanhADAA::~TanhADAA() {
  RTFree(mWorld, osBuffer);
  RTFree(mWorld, cpyBuf);
}

inline double TanhADAA::tanh(const double &v) { return std::tanh(v); }

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

  const int osBufferSize = nSamples * fScale;
  // println("Starting a next_aa call");
  const float *sig = in(Input);
  float *outbuf = out(Out1);
  const int adlevel = static_cast<int>(in0(AntiDerivativeLevel));

  // up sample -- results are stored in osBuffer
  os.processSamplesUp(sig, cpyBuf, osBuffer);

  // process
  for (auto i = 0; i < osBufferSize; ++i) {
    if (adlevel == ADAA::AntiDerivativeLevel::First) {
      double curr = ADAA::next_first_adaa(
          osBuffer[i], x1, ad1_x1, TanhADAA::tanh, TanhADAA::tanh_first_ad);
      osBuffer[i] = curr;

    } else if (adlevel == ADAA::AntiDerivativeLevel::Second) {
      double curr = ADAA::next_second_adaa(
          osBuffer[i], x1, x2, ad2_x0, ad2_x1, d2, TanhADAA::tanh,
          TanhADAA::tanh_first_ad, TanhADAA::tanh_second_ad);
      osBuffer[i] = curr;
    }
  }

  // down sample -- results stored in outbuf
  os.processSamplesDown(outbuf, osBuffer);
}

}  // namespace JSCDSP::TanhADAA

PluginLoad(HardClipADAAUGens) {
  // Plugin magic
  ft = inTable;
  registerUnit<JSCDSP::HardClipADAA::HardClipADAA>(ft, "HardClipADAA", false);
  registerUnit<JSCDSP::TanhADAA::TanhADAA>(ft, "TanhADAA", false);
}
