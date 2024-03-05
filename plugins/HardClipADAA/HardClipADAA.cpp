// Nonlinear shapers with Anti Derivative Anti Aliasing 
// James Squires (squires.jr@gmail.com)
//
// Adapted from Chow DSP

#include "SC_PlugIn.hpp"
#include "HardClipADAA.hpp"
#include "SC_PlugIn.h"
#include <algorithm>
#include <cmath>
#include "Utils.hpp"
#include "Li2.hpp"


static InterfaceTable* ft;

namespace TOL {
  constexpr double TOL = 1.0e-5;
}

namespace ADAA {

template<typename Func, typename Ad_Func> 
static inline double next_first_adaa(const double& s, double& x1, double& ad1_x1, Func f, Ad_Func first_ad) 
{
  double res;
  double ad1_x = first_ad(s);
  double diff = s - x1;

  if (std::abs(diff) <= TOL::TOL) {
    res = f(0.5 * (s + x1));
  } else {
    res = (ad1_x - ad1_x1) / diff;
  }

  ad1_x1 = ad1_x;
  x1 = s;

  return res;
}

template<typename Func, typename FuncFirstAD, typename FuncSecondAD>
inline double next_second_adaa(const double& s, 
                                      double& x1, double& x2, 
                                      double& ad2_x0, double& ad2_x1, 
                                      double& d2, 
                                      Func f, 
                                      FuncFirstAD f_first_ad, FuncSecondAD f_second_ad) 
{
  double res;
  double d1 = calcD(s, x1, ad2_x0, ad2_x1, f_first_ad, f_second_ad);

  // x_n - x_{n-1} <= epsilon, then use 
  if (std::abs(s - x2) <= TOL::TOL) { // fallback
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


template<typename FuncFirstAD, typename FuncSecondAD>
inline double calcD(const double& v0,
                                         double& x1, 
                                         double& ad2_x0, double& ad2_x1,
                                         FuncFirstAD f_first_ad, FuncSecondAD f_second_ad) 
{
  ad2_x0 = f_second_ad(v0);

  if (std::abs(v0 - x1) <= TOL::TOL) {
    return f_first_ad(0.5 * (v0 + x1));
  } else {
    return (ad2_x0 - ad2_x1) / (v0 - x1);
  }

}


template<typename Func, typename FuncFirstAD, typename FuncSecondAD>
inline double fallback(const double& x,
                              double& x1, double& x2,
                              double& ad2_x1,
                              Func f, FuncFirstAD f_first_ad, FuncSecondAD f_second_ad) 
{
  const double xbar = 0.5 * (x + x2);
  const double delta = xbar - x1;

  if (std::abs(delta) <= TOL::TOL) {
    return f(0.5 * (xbar + x1));
  } else {
    return (2.0 / delta) * (f_first_ad(xbar) + ((ad2_x1 - f_second_ad(xbar)) / delta));
  }
}


} // Namespace ADAA




namespace HardClipADAA {

HardClipADAA::HardClipADAA() {


  os.init(2, (int)inBufferSize(0));

  fcoefs.gen_coefs((int)sampleRate());

  acoeff = fcoefs.acoeff;
  bcoeff = fcoefs.bcoeff;
  gain = *fcoefs.gain;

  xv1 = (double*)RTAlloc(mWorld, (npole + 1) * sizeof(double));
  yv1 = (double*)RTAlloc(mWorld, (npole + 1) * sizeof(double));
  xv2 = (double*)RTAlloc(mWorld, (npole + 1) * sizeof(double));
  yv2 = (double*)RTAlloc(mWorld, (npole + 1) * sizeof(double));

  osBuffer = (double*)RTAlloc(mWorld, inBufferSize(0) * 2 * sizeof(double));

  assert(xv1 != NULL);
  assert(yv1 != NULL);
  assert(xv2 != NULL);
  assert(yv2 != NULL);
  assert(osBuffer != NULL);

  std::fill(&xv1[0], &xv1[npole + 1], 0.0);
  std::fill(&yv1[0], &yv1[npole + 1], 0.0);
  std::fill(&xv2[0], &xv2[npole + 1], 0.0);
  std::fill(&yv2[0], &yv2[npole + 1], 0.0);

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
  RTFree(mWorld, xv1);
  RTFree(mWorld, yv1);
  RTFree(mWorld, xv2);
  RTFree(mWorld, yv2);
  RTFree(mWorld, osBuffer);
}


inline double HardClipADAA::filter_next(double v, double* xv, double* yv) {
  double out = 0.0;
  int j;

  for (j = 0; j < nzero; j++) {
    xv[j] = xv[j+1];
  }
  xv[nzero] = v / gain;

  for (j = 0; j < npole; j++) {
    yv[j] = yv[j + 1];
  }

  for (j = 0; j <= nzero; j++) {
    out += xv[j] * bcoeff[j];
  }

  for (j = 0; j < npole; j++) {
    out -= yv[j] * acoeff[j];
  }
  yv[npole] = out;

  return out;

}

inline double HardClipADAA::clip(const double& v) {
  return (std::abs(v) <= 1.0f) ? v : 1.0;
}

inline double HardClipADAA::hc_first_ad(const double& v) {
  return (std::abs(v) <= 1.0f) ? (v * v) * 0.5 : (v * signum(v)) - 0.5;
}

inline double HardClipADAA::hc_second_ad(const double& v) {
  if (std::abs(v) <= 1.0f) {
    return std::pow(v, 3) / 6.0;
  } else {
    return (((std::pow(v, 2) * 0.5) + (1.0/6.0)) * signum(v)) - (v * 0.5);
  }
}

void HardClipADAA::next_aa(int nSamples) {
  int i;
  int twoNSamples = nSamples * 2;

  const float* sig = in(Input);

  float* outbuf = out(Out1); const int adlevel = (int)in0(AntiDerivativeLevel);

  // reset Oversampling for oversampling  
  os.reset();
  // upsample
  os.processSamplesUp(sig, nSamples);

  // process
  for (double s : os.getProcessedSamples()) {
    if (adlevel == ADAA::AntiDerivativeLevel::FirstOrder){
      osBuffer[i] = ADAA::next_first_adaa(s, 
                                          x1, 
                                          ad1_x1, 
                                          HardClipADAA::clip, 
                                          HardClipADAA::hc_first_ad);

    } else if (adlevel == ADAA::AntiDerivativeLevel::SecondOrder) {
      osBuffer[i] = ADAA::next_second_adaa(s,
                                           x1, x2,
                                           ad2_x0, ad2_x1, d2,
                                           HardClipADAA::clip,
                                           HardClipADAA::hc_first_ad,
                                           HardClipADAA::hc_second_ad);

    }
  }
  
  // filter again
  os.processSamplesDown();

  for (double s : os.getProcessedSamples()) {
    outbuf[i] = static_cast<float>(s);
  }

}

} // Namespace HardClipADAA




namespace TanhADAA {

TanhADAA::TanhADAA() {

  fcoefs.gen_coefs((int)sampleRate());

  acoeff = fcoefs.acoeff;
  bcoeff = fcoefs.bcoeff;
  gain = *fcoefs.gain;

  xv1 = (double*)RTAlloc(mWorld, (npole + 1) * sizeof(double));
  yv1 = (double*)RTAlloc(mWorld, (npole + 1) * sizeof(double));
  xv2 = (double*)RTAlloc(mWorld, (npole + 1) * sizeof(double));
  yv2 = (double*)RTAlloc(mWorld, (npole + 1) * sizeof(double));

  osBuffer = (double*)RTAlloc(mWorld, inBufferSize(0) * 2 * sizeof(double));

  assert(xv1 != NULL);
  assert(yv1 != NULL);
  assert(xv2 != NULL);
  assert(yv2 != NULL);
  assert(osBuffer != NULL);

  std::fill(&xv1[0], &xv1[npole + 1], 0.0);
  std::fill(&yv1[0], &yv1[npole + 1], 0.0);
  std::fill(&xv2[0], &xv2[npole + 1], 0.0);
  std::fill(&yv2[0], &yv2[npole + 1], 0.0);

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
  RTFree(mWorld, xv1);
  RTFree(mWorld, yv1);
  RTFree(mWorld, xv2);
  RTFree(mWorld, yv2);
  RTFree(mWorld, osBuffer);
}


inline double TanhADAA::filter_next(double v, double* xv, double* yv) {
  double out = 0.0;
  int j;

  for (j = 0; j < nzero; j++) {
    xv[j] = xv[j+1];
  }
  xv[nzero] = v / gain;

  for (j = 0; j < npole; j++) {
    yv[j] = yv[j + 1];
  }

  for (j = 0; j <= nzero; j++) {
    out += xv[j] * bcoeff[j];
  }

  for (j = 0; j < npole; j++) {
    out -= yv[j] * acoeff[j];
  }
  yv[npole] = out;

  return out;

}

inline double TanhADAA::tanh_first_ad(const double& v) {
  return std::log (std::cosh (v));
}

inline double TanhADAA::tanh_second_ad(const double& v) {
  const auto exp = std::exp (-2 * v);
  return 0.5 * ((double) polylogarithm::Li2 (-exp) 
          - v * (v + 2.0 * std::log (exp + 1.) - 2.0 * std::log (std::cosh (v))))
          + (Constants::PI_SQRD / 24.0);
}

void TanhADAA::next_aa(int nSamples) {
  int i;
  int twoNSamples = nSamples * 2;

  const float* sig = in(Input);

  float* outbuf = out(Out1); const int adlevel = (int)in0(AntiDerivativeLevel);
  // upsample 
  for (i = 0; i < nSamples; ++i) {
    osBuffer[i << 1] = static_cast<double>(sig[i]);
    osBuffer[(i << 1) + 1] = 0.0;
  }

  // filter
  for (i = 0; i < twoNSamples; ++i) {
    osBuffer[i] = filter_next(osBuffer[i], xv1, yv1);
  }

  // process
  for (i = 0; i < twoNSamples; ++i) {

    if (adlevel == ADAA::AntiDerivativeLevel::FirstOrder){

      osBuffer[i] = ADAA::next_first_adaa(osBuffer[i], 
                                          x1, 
                                          ad1_x1, 
                                          [](double v) { return std::tanh(v); }, 
                                          TanhADAA::tanh_first_ad);

    } else if (adlevel == ADAA::AntiDerivativeLevel::SecondOrder) {

      osBuffer[i] = ADAA::next_second_adaa(osBuffer[i],
                                           x1, x2,
                                           ad2_x0, ad2_x1, d2,
                                           [](double v) { return std::tanh(v); },
                                           TanhADAA::tanh_first_ad,
                                           TanhADAA::tanh_second_ad);

    }

  }

  // filter again
  for (i = 0; i < twoNSamples; ++i) {
    osBuffer[i] = filter_next(osBuffer[i], xv2, yv2);
  }

  // downsample and out
  for (i = 0; i < nSamples; ++i) {
    outbuf[i] = static_cast<float>(osBuffer[i << 1]);
  }


}

} // namespace TanhADAA



PluginLoad(HardClipADAAUGens) {
  // Plugin magic
  ft = inTable;
  registerUnit<HardClipADAA::HardClipADAA>(ft, "HardClipADAA", false);
  registerUnit<TanhADAA::TanhADAA>(ft, "TanhADAA", false);
}
