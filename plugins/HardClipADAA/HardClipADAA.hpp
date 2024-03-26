// Copyright 2024 James Squires
// PluginHardClipADAA.hpp
// James Squires (squires.jr@gmail.com)

#pragma once

#include "Oversampling.hpp"
#include "SC_PlugIn.hpp"

namespace JSCDSP::ADAA {

template <typename Func, typename Ad_Func>
static inline double next_first_adaa(const double& s, double* x1,
                                     double* ad1_x1, Func f, Ad_Func first_ad);

template <typename Func, typename FuncFirstAD, typename FuncSecondAD>
inline double next_second_adaa(const double& s, double* x1, double* x2,
                               double* ad2_x0, double* ad2_x1, double* d2,
                               Func f, FuncFirstAD f_first_ad,
                               FuncSecondAD f_second_ad);

template <typename FuncFirstAD, typename FuncSecondAD>
inline double calcD(const double& v0, const double& x1, double* ad2_x0,
                    const double& ad2_x1, FuncFirstAD f_first_ad,
                    FuncSecondAD f_second_ad);

template <typename Func, typename FuncFirstAD, typename FuncSecondAD>
inline double fallback(const double& x, const double& x1, const double& x2,
                       const double& ad2_x1, Func f, FuncFirstAD f_first_ad,
                       FuncSecondAD f_second_ad);

}  // namespace JSCDSP::ADAA

namespace JSCDSP::HardClipADAA {

class HardClipADAA : public SCUnit {
 public:
  HardClipADAA();
  // Destructor
  ~HardClipADAA();

 private:
  // Calc function
  void next_aa(int nSamples);
  static inline double signum(const double& v);
  static inline double clip(const double& v);
  static inline double hc_first_ad(const double& v);
  static inline double hc_second_ad(const double& v);

  double* osBuffer;
  double* cpyBuf;
  Oversampling::Oversampling os;


  int fScale{2};
  double x1{0.0f};
  double ad1_x1{0.0f};
  double x2{0.0f};
  double d2{0.0f};
  double ad2_x1{0.0f};
  double ad2_x0{0.0f};

  enum InputParams { Input, AntiDerivativeLevel, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };
};

}  // namespace JSCDSP::HardClipADAA

namespace JSCDSP::TanhADAA {

class TanhADAA : public SCUnit {
 public:
  TanhADAA();
  // Destructor
  ~TanhADAA();

 private:
  // Calc function
  void next_aa(int nSamples);
  static inline double tanh(const double& v);
  static inline double tanh_first_ad(const double& v);
  static inline double tanh_second_ad(const double& v);

  double* osBuffer;
  double* cpyBuf;

  int fScale{2};
  Oversampling::Oversampling os;

  double x1{0.0};
  double ad1_x1{0.0};
  double x2{0.0};
  double d2{0.0};
  double ad2_x1{0.0};
  double ad2_x0{0.0};

  enum InputParams { Input, AntiDerivativeLevel, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };
};

}  // namespace JSCDSP::TanhADAA
