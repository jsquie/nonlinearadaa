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

constexpr int FILTER_TAP_NUM = 209;

static const double filter_taps[FILTER_TAP_NUM] = {
    -0.0016550762341430308,  -0.0041665939310652385,  -0.003594191447486389,
    0.00251719066089494,     0.009204232105743072,    0.00891685978767867,
    0.0025929990340545586,   -0.0009241688686624757,  0.001780417290846837,
    0.0040583237055417075,   0.0012718880881058892,   -0.0015563530079133385,
    0.0002721535795504213,   0.0025035930967126315,   0.0005429619868192588,
    -0.0017453123023537607,  0.000022984138603901914, 0.0020537985664696805,
    0.00013860481511378278,  -0.0018438080070986866,  0.00013467879098560712,
    0.001974941764859404,    -0.00018853242904265353, -0.0019542465345495637,
    0.0003762055310051086,   0.0020296451360248486,   -0.0005169276482732894,
    -0.0020684792621014857,  0.0006994817682676137,   0.0021265551428543772,
    -0.0008822257189953905,  -0.0021674203444939926,  0.0010935813680486932,
    0.0022115890943483034,   -0.001314750098015205,   -0.0022412051839012007,
    0.0015560490952100908,   0.0022569517475881764,   -0.0018209974880386698,
    -0.002266301757028872,   0.0020937340035129225,   0.0022477597114187653,
    -0.0023927521827856475,  -0.0022146568586112986,  0.0027050436381634005,
    0.002156314932451038,    -0.003031065300689767,   -0.002064650858536239,
    0.0033792643077078716,   0.0019474808728353638,   -0.0037376568321684455,
    -0.0017910129887243696,  0.004112825480457824,    0.0015960522407615228,
    -0.004502958711563328,   -0.0013612537988873558,  0.004903329017739637,
    0.001079361495260796,    -0.005313440162938722,   -0.000745420496662418,
    0.00573623700664673,     0.00036334614253456674,  -0.006156309094600619,
    0.00008808856155454664,  0.006589753951165301,    -0.0005956547623096374,
    -0.007014724209099923,   0.0011845526388683388,   0.007445042007534717,
    -0.0018497468084161031,  -0.007868705069367566,   0.0026066792019242466,
    0.008285357331304153,    -0.0034688911024595,     -0.008698582991922315,
    0.0044411237414403925,   0.009092986565733343,    -0.005557414337636592,
    -0.009483553419766278,   0.006823202063313745,    0.009848519291251277,
    -0.008289679120818303,   -0.010201859443027823,   0.009988471477688032,
    0.01053012360277123,     -0.011992197126080751,   -0.010834696922933549,
    0.014396832948754698,    0.011115959799250405,    -0.017345990204051556,
    -0.01136588090176962,    0.021087101188298568,    0.01158952167903816,
    -0.026030444783605375,   -0.011780011049371056,   0.03296919343261422,
    0.011939078466010828,    -0.04359427254454018,    -0.012063942524891857,
    0.0623135104924757,      0.012154411619644348,    -0.10529117509190798,
    -0.01220839945864449,    0.3180391210647029,      0.5122272478248557,
    0.3180391210647029,      -0.01220839945864449,    -0.10529117509190798,
    0.012154411619644348,    0.0623135104924757,      -0.012063942524891857,
    -0.04359427254454018,    0.011939078466010828,    0.03296919343261422,
    -0.011780011049371056,   -0.026030444783605375,   0.01158952167903816,
    0.021087101188298568,    -0.01136588090176962,    -0.017345990204051556,
    0.011115959799250405,    0.014396832948754698,    -0.010834696922933549,
    -0.011992197126080751,   0.01053012360277123,     0.009988471477688032,
    -0.010201859443027823,   -0.008289679120818303,   0.009848519291251277,
    0.006823202063313745,    -0.009483553419766278,   -0.005557414337636592,
    0.009092986565733343,    0.0044411237414403925,   -0.008698582991922315,
    -0.0034688911024595,     0.008285357331304153,    0.0026066792019242466,
    -0.007868705069367566,   -0.0018497468084161031,  0.007445042007534717,
    0.0011845526388683388,   -0.007014724209099923,   -0.0005956547623096374,
    0.006589753951165301,    0.00008808856155454664,  -0.006156309094600619,
    0.00036334614253456674,  0.00573623700664673,     -0.000745420496662418,
    -0.005313440162938722,   0.001079361495260796,    0.004903329017739637,
    -0.0013612537988873558,  -0.004502958711563328,   0.0015960522407615228,
    0.004112825480457824,    -0.0017910129887243696,  -0.0037376568321684455,
    0.0019474808728353638,   0.0033792643077078716,   -0.002064650858536239,
    -0.003031065300689767,   0.002156314932451038,    0.0027050436381634005,
    -0.0022146568586112986,  -0.0023927521827856475,  0.0022477597114187653,
    0.0020937340035129225,   -0.002266301757028872,   -0.0018209974880386698,
    0.0022569517475881764,   0.0015560490952100908,   -0.0022412051839012007,
    -0.001314750098015205,   0.0022115890943483034,   0.0010935813680486932,
    -0.0021674203444939926,  -0.0008822257189953905,  0.0021265551428543772,
    0.0006994817682676137,   -0.0020684792621014857,  -0.0005169276482732894,
    0.0020296451360248486,   0.0003762055310051086,   -0.0019542465345495637,
    -0.00018853242904265353, 0.001974941764859404,    0.00013467879098560712,
    -0.0018438080070986866,  0.00013860481511378278,  0.0020537985664696805,
    0.000022984138603901914, -0.0017453123023537607,  0.0005429619868192588,
    0.0025035930967126315,   0.0002721535795504213,   -0.0015563530079133385,
    0.0012718880881058892,   0.0040583237055417075,   0.001780417290846837,
    -0.0009241688686624757,  0.0025929990340545586,   0.00891685978767867,
    0.009204232105743072,    0.00251719066089494,     -0.003594191447486389,
    -0.0041665939310652385,  -0.0016550762341430308};


  // M = FIRFilter::FILTER_TAP_NUM;
  M = 209;

  // assert (filter_taps != nullptr);

  osBuffer = reinterpret_cast<double *>(
      RTAlloc(mWorld, inBufferSize(0) * in0(OverSample) * sizeof(double)));

  float Mf = static_cast<float>(M);
  int first_buffer_size = std::ceil(Mf / in0(OverSample));
  int last_buffer_size = std::floor(Mf / in0(OverSample));

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
