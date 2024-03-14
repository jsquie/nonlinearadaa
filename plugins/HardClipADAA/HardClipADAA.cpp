// Copyright 2024 James Squires
// Nonlinear shapers with Anti Derivative Anti Aliasing
// James Squires (squires.jr@gmail.com)
//
// Adapted from Chow DSP
#include <cmath>
#include <iostream>

#include "FIRFilter.hpp"
#include "HardClipADAA.hpp"
#include "Li2.hpp"
#include "SC_InterfaceTable.h"
#include "SC_PlugIn.hpp"
#include "Utils.hpp"

static InterfaceTable *ft;

static unsigned int numCalls = 0;

namespace JSCDSP::ADAA {

constexpr double TOL = 1.0e-5;

template <typename Func, typename Ad_Func>
static inline double next_first_adaa(const double &s, double *x1,
                                     double *ad1_x1, Func f, Ad_Func first_ad) {
  double res;
  double ad1_x = first_ad(s);
  double diff = s - *x1;

  if (std::abs(diff) <= TOL) {
    res = f(0.5 * (s + *x1));
  } else {
    res = (ad1_x - *ad1_x1) / diff;
  }

  *ad1_x1 = ad1_x;
  *x1 = s;

  return res;
}

template <typename Func, typename FuncFirstAD, typename FuncSecondAD>
inline double next_second_adaa(const double &s, double *x1, double *x2,
                               double *ad2_x0, double *ad2_x1, double *d2,
                               Func f, FuncFirstAD f_first_ad,
                               FuncSecondAD f_second_ad) {
  double res;
  double d1 = calcD(s, *x1, ad2_x0, *ad2_x1, f_first_ad, f_second_ad);

  // x_n - x_{n-1} <= epsilon, then use
  if (std::abs(s - *x2) <= TOL) {  // fallback
    res = fallback(s, *x1, *x2, *ad2_x1, f, f_first_ad, f_second_ad);
  } else {
    res = (2.0 / (s - *x2)) * (d1 - *d2);
  }

  *d2 = d1;
  *x2 = *x1;
  *x1 = s;
  *ad2_x1 = *ad2_x0;

  return res;
}

template <typename FuncFirstAD, typename FuncSecondAD>
inline double calcD(const double &v0, const double &x1, double *ad2_x0,
                    const double &ad2_x1, FuncFirstAD f_first_ad,
                    FuncSecondAD f_second_ad) {
  *ad2_x0 = f_second_ad(v0);

  if (std::abs(v0 - x1) <= TOL) {
    return f_first_ad(0.5 * (v0 + x1));
  } else {
    return (*ad2_x0 - ad2_x1) / (v0 - x1);
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

enum AntiDerivativeLevel { First = 1, Second = 2 };

}  // namespace JSCDSP::ADAA

namespace JSCDSP::HardClipADAA {

HardClipADAA::HardClipADAA() {
  // initialize with factor oversample and buffersize
  // define stopband atten
  double alpha = 20.0;
  // define transition bandwidth between pass and stop bands
  double transition_bw = 0.1;
  // determine length of kernel -- M
  M = FIRFilter::computeKaiserLength<double>(transition_bw, alpha);
  double beta = FIRFilter::computeKaiserBeta<double>(alpha);

  M = (M % 2 == 0) ? ++M : M;
  M = 329;

  std::cout << "M: " << M << std::endl;

  // allocate resources for os buffer
  osBuffer = reinterpret_cast<double *>(
      RTAlloc(mWorld, inBufferSize(0) * in0(OverSample) * sizeof(double)));

  os.setOsBuffer(osBuffer);
  // allocate resources for kaiser window (size double * M)
  // allocate resources for filter kernel (size double * M)
  // allocate resrouces for circ buffer (size double * M)
  // allocate resources for filter output ( ???? )
  kWindowBuf = reinterpret_cast<double *>(RTAlloc(mWorld, M * sizeof(double)));

  fKernelBuf = reinterpret_cast<double *>(RTAlloc(mWorld, M * sizeof(double)));

  std::cout << "Resources are allocated. Proceeding to design filter"
            << std::endl;
  // design filter -- results in kernel stored in fKernelBuf
  /*

  FIR filter designed with
  http://t-filter.appspot.com

  sampling frequency: 88200 Hz

  * 0 Hz - 22000 Hz
    gain = 1
    desired ripple = 1 dB
    actual ripple = 0.7405780597190046 dB

  * 23000 Hz - 44100 Hz
    gain = 0
    desired attenuation = -120 dB
    actual attenuation = -120.13314207825391 dB

  */

#define FILTER_TAP_NUM 329
  M = FILTER_TAP_NUM;
  static double filter_taps[FILTER_TAP_NUM] = {
      0.00005955536388890257,  0.0005086273557871056,   0.002002845673048492,
      0.004819098163541085,    0.007633910679898034,    0.007665506879457962,
      0.003465998685522977,    -0.0022382534565408567,  -0.0043163738936135105,
      -0.0011524093722430143,  0.0027672788957769077,   0.0023529563038693305,
      -0.0012392877502301698,  -0.002490682346558991,   0.00016591921970796482,
      0.0022199111232029363,   0.0005046446566548329,   -0.0018352400372817139,
      -0.000890593251161442,   0.0014603875326106763,   0.0010920931586451198,
      -0.0011440724158604866,  -0.0011875218507201156,  0.0008924668056188886,
      0.0012289068113610897,   -0.0006890201021156335,  -0.0012338249085039408,
      0.0005290357996686363,   0.001214377051757208,    -0.0004154808926984674,
      -0.0011933400723739818,  0.00033408095066981164,  0.001179723712785923,
      -0.0002674962473332636,  -0.0011631801500141643,  0.0002216672511884556,
      0.0011502744366912003,   -0.00019532929152349252, -0.0011515887833740614,
      0.00017026628738050264,  0.0011499921157861709,   -0.0001588078721630991,
      -0.0011517489739520548,  0.00016489912554385066,  0.0011740065342130393,
      -0.00016277759510698479, -0.0011920082005782149,  0.00016925516782108883,
      0.0012128735599406444,   -0.00018658237571223123, -0.0012480319480878664,
      0.00019468261742233773,  0.0012699152863328684,   -0.0002272558739618597,
      -0.0013163830136113173,  0.0002454416022234823,   0.0013509120275368252,
      -0.00028176122778699944, -0.0014023547424380103,  0.00031041389101307213,
      0.0014487801055858285,   -0.0003475187661534844,  -0.0014995488164712089,
      0.00038885658924687393,  0.0015522532185385596,   -0.00043772728610685946,
      -0.0016117673018712854,  0.0004886482593258717,   0.0016735268648443598,
      -0.0005442311650481099,  -0.0017384539494585789,  0.0006040607324215597,
      0.0018062372879560987,   -0.0006672561408347555,  -0.0018741590017286515,
      0.0007366258181371059,   0.0019423975489675505,   -0.0008152808109399897,
      -0.0020156502802753545,  0.0008987017547054996,   0.002090309533453088,
      -0.0009890634948111381,  -0.0021667294459222094,  0.0010873448957003631,
      0.002246512906620246,    -0.001191024370009242,   -0.002325829676405825,
      0.0013040369692723115,   0.0024069715880707227,   -0.0014256601695837423,
      -0.002489534568009997,   0.001555790147669491,    0.002571514990956858,
      -0.001698278667799383,   -0.0026570496676622973,  0.0018493497907988994,
      0.0027423201067591726,   -0.00201215218634822,    -0.0028282298964724764,
      0.0021881261695355213,   0.002917661506350236,    -0.002372414103372207,
      -0.0030028827653562594,  0.002573752130853033,    0.0030905799441261635,
      -0.002787092574359089,   -0.0031743503895247504,  0.0030202176423288415,
      0.0032602543572685513,   -0.0032692763505512224,  -0.0033437051203688306,
      0.0035402840987689634,   0.0034294860433242455,   -0.0038301208859661774,
      -0.0035136548516605534,  0.004142903803499536,    0.003596841210119382,
      -0.0044805095776605545,  -0.003678110561878164,   0.004846504817964823,
      0.00375766110973799,     -0.0052440538322135405,  -0.003834777863993236,
      0.005678387865873917,    0.003911158179364791,    -0.006152023204316004,
      -0.003984066411633326,   0.006673468868571328,    0.004055889106924824,
      -0.007247577326601133,   -0.00412447118007218,    0.007884547568392372,
      0.004189783464988519,    -0.008596500332184301,   -0.004252966998921284,
      0.009396328509984515,    0.0043120950018453634,   -0.010304334208457321,
      -0.004367589828556944,   0.011345767765827013,    0.004419981225170736,
      -0.01255360999257909,    -0.004468266595084712,   0.013975602922147494,
      0.004513484006735141,    -0.015677221551188242,   -0.004555097322056134,
      0.017756198025810636,    0.0045929017019218496,   -0.020362571712357058,
      -0.004627164215982237,   0.023737718550898446,    0.0046562662701334955,
      -0.02830214600831298,    -0.004682016258796661,   0.03484465090244982,
      0.004701178931061656,    -0.04506446236482974,    -0.004716923404808415,
      0.06336794648719264,     0.004726159133302609,    -0.10592639655710469,
      -0.0047319810673935574,  0.31825072417452005,     0.5047337174872589,
      0.31825072417452005,     -0.0047319810673935574,  -0.10592639655710469,
      0.004726159133302609,    0.06336794648719264,     -0.004716923404808415,
      -0.04506446236482974,    0.004701178931061656,    0.03484465090244982,
      -0.004682016258796661,   -0.02830214600831298,    0.0046562662701334955,
      0.023737718550898446,    -0.004627164215982237,   -0.020362571712357058,
      0.0045929017019218496,   0.017756198025810636,    -0.004555097322056134,
      -0.015677221551188242,   0.004513484006735141,    0.013975602922147494,
      -0.004468266595084712,   -0.01255360999257909,    0.004419981225170736,
      0.011345767765827013,    -0.004367589828556944,   -0.010304334208457321,
      0.0043120950018453634,   0.009396328509984515,    -0.004252966998921284,
      -0.008596500332184301,   0.004189783464988519,    0.007884547568392372,
      -0.00412447118007218,    -0.007247577326601133,   0.004055889106924824,
      0.006673468868571328,    -0.003984066411633326,   -0.006152023204316004,
      0.003911158179364791,    0.005678387865873917,    -0.003834777863993236,
      -0.0052440538322135405,  0.00375766110973799,     0.004846504817964823,
      -0.003678110561878164,   -0.0044805095776605545,  0.003596841210119382,
      0.004142903803499536,    -0.0035136548516605534,  -0.0038301208859661774,
      0.0034294860433242455,   0.0035402840987689634,   -0.0033437051203688306,
      -0.0032692763505512224,  0.0032602543572685513,   0.0030202176423288415,
      -0.0031743503895247504,  -0.002787092574359089,   0.0030905799441261635,
      0.002573752130853033,    -0.0030028827653562594,  -0.002372414103372207,
      0.002917661506350236,    0.0021881261695355213,   -0.0028282298964724764,
      -0.00201215218634822,    0.0027423201067591726,   0.0018493497907988994,
      -0.0026570496676622973,  -0.001698278667799383,   0.002571514990956858,
      0.001555790147669491,    -0.002489534568009997,   -0.0014256601695837423,
      0.0024069715880707227,   0.0013040369692723115,   -0.002325829676405825,
      -0.001191024370009242,   0.002246512906620246,    0.0010873448957003631,
      -0.0021667294459222094,  -0.0009890634948111381,  0.002090309533453088,
      0.0008987017547054996,   -0.0020156502802753545,  -0.0008152808109399897,
      0.0019423975489675505,   0.0007366258181371059,   -0.0018741590017286515,
      -0.0006672561408347555,  0.0018062372879560987,   0.0006040607324215597,
      -0.0017384539494585789,  -0.0005442311650481099,  0.0016735268648443598,
      0.0004886482593258717,   -0.0016117673018712854,  -0.00043772728610685946,
      0.0015522532185385596,   0.00038885658924687393,  -0.0014995488164712089,
      -0.0003475187661534844,  0.0014487801055858285,   0.00031041389101307213,
      -0.0014023547424380103,  -0.00028176122778699944, 0.0013509120275368252,
      0.0002454416022234823,   -0.0013163830136113173,  -0.0002272558739618597,
      0.0012699152863328684,   0.00019468261742233773,  -0.0012480319480878664,
      -0.00018658237571223123, 0.0012128735599406444,   0.00016925516782108883,
      -0.0011920082005782149,  -0.00016277759510698479, 0.0011740065342130393,
      0.00016489912554385066,  -0.0011517489739520548,  -0.0001588078721630991,
      0.0011499921157861709,   0.00017026628738050264,  -0.0011515887833740614,
      -0.00019532929152349252, 0.0011502744366912003,   0.0002216672511884556,
      -0.0011631801500141643,  -0.0002674962473332636,  0.001179723712785923,
      0.00033408095066981164,  -0.0011933400723739818,  -0.0004154808926984674,
      0.001214377051757208,    0.0005290357996686363,   -0.0012338249085039408,
      -0.0006890201021156335,  0.0012289068113610897,   0.0008924668056188886,
      -0.0011875218507201156,  -0.0011440724158604866,  0.0010920931586451198,
      0.0014603875326106763,   -0.000890593251161442,   -0.0018352400372817139,
      0.0005046446566548329,   0.0022199111232029363,   0.00016591921970796482,
      -0.002490682346558991,   -0.0012392877502301698,  0.0023529563038693305,
      0.0027672788957769077,   -0.0011524093722430143,  -0.0043163738936135105,
      -0.0022382534565408567,  0.003465998685522977,    0.007665506879457962,
      0.007633910679898034,    0.004819098163541085,    0.002002845673048492,
      0.0005086273557871056,   0.00005955536388890257};
  for (int i = 0; i < M; ++i) {
    fKernelBuf[i] = filter_taps[i];
  }

  os.init(in0(OverSample), inBufferSize(0), M, fKernelBuf);

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
  RTFree(mWorld, kWindowBuf);
  RTFree(mWorld, fKernelBuf);
}

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

  // filter with FIR filter
  // FIRFilter::linear_convolve(sig, fKernelBuf, *cBuf, fOutput, nSamples, M);

  double y0 = 0.0f;
  double y02 = 0.0f;
  double y1 = 0.0f;
  double y12 = 0.0f;

  os.processSamplesUp(sig);

  // process
  for (auto s : os.getProcessedSamples()) {
    if (adlevel == JSCDSP::ADAA::AntiDerivativeLevel::First) {
      s = JSCDSP::ADAA::next_first_adaa(s, &x1, &ad1_x1,
                                                  HardClipADAA::clip,
                                                  HardClipADAA::hc_first_ad);

    } else if (adlevel == ADAA::AntiDerivativeLevel::Second) {
      s = JSCDSP::ADAA::next_second_adaa(
          s, &x1, &x2, &ad2_x0, &ad2_x1, &d2, HardClipADAA::clip,
          HardClipADAA::hc_first_ad, HardClipADAA::hc_second_ad);
    }
  }

  os.processSamplesDown(outbuf);
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
