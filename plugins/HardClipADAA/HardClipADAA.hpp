// PluginHardClipADAA.hpp
// James Squires (squires.jr@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"
#include "IIRFilterCoefs.hpp"
#include "Oversampling.hpp"

namespace ADAA {

enum AntiDerivativeLevel { FirstOrder = 1, SecondOrder = 2 };

template<typename Func, typename Ad_Func> 
static inline double next_first_adaa(const double& s, 
                                     double*& x1, double*& ad1_x1, 
                                     Func f, Ad_Func first_ad);

template<typename Func, typename FuncFirstAD, typename FuncSecondAD>
inline double next_second_adaa(const double& s, 
                              double*& x1, double*& x2, 
                              double*& ad2_x0, double*& ad2_x1, 
                              double*& d2, 
                              Func f, 
                              FuncFirstAD f_first_ad, FuncSecondAD f_second_ad);

template<typename FuncFirstAD, typename FuncSecondAD>
inline double calcD(const double& v0,
                   double*& x1, 
                   double*& ad2_x0, double*& ad2_x1,
                   FuncFirstAD f_first_ad, FuncSecondAD f_second_ad); 

template<typename Func, typename FuncFirstAD, typename FuncSecondAD>
inline double fallback(const double& x,
                      double*& x1, double*& x2,
                      double*& ad2_x1,
                      Func f, FuncFirstAD f_first_ad, FuncSecondAD f_second_ad);



} // Namespace ADAA




namespace HardClipADAA {

class HardClipADAA : public SCUnit {

public:
  HardClipADAA();
  // Destructor
  ~HardClipADAA();
  IIRFilterCoefs fcoefs;
  // Oversampling::Oversampling os;

private:
  // Calc function
  static inline double signum(const double& v) { return v < 0 ? -1.0 : 1.0; };
  void next_aa(int nSamples);
  static inline double clip(const double& v);
  static inline double hc_first_ad(const double& v);
  static inline double hc_second_ad(const double& v);
  inline double filter_next(double v, double* xv, double* yv);

  double x1;
  double ad1_x1;
  double x2;
  double d2;
  double ad2_x1;
  double ad2_x0;

  double* osBuffer;

  const int npole = 4;
  const int nzero = 4;
  // for samplerate == 48000
  const double* acoeff;
  const double* bcoeff;
  double gain;

  double* xv1;
  double* yv1;
  double* xv2;
  double* yv2;

  float* dbup;
  float* down1;
  float* dbdown;

  Oversampling::Oversampling os;

  enum InputParams { Input, AntiDerivativeLevel, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };


};

} // namespace HardClipADAA

namespace TanhADAA {
 
class TanhADAA : public SCUnit {

public:
  TanhADAA();
  // Destructor
  ~TanhADAA();
  IIRFilterCoefs fcoefs;

private:
  // Calc function
  void next_aa(int nSamples);
  static inline double tanh_first_ad(const double& v);
  static inline double tanh_second_ad(const double& v);
  inline double filter_next(double v, double* xv, double* yv);

  double x1;
  double ad1_x1;
  double x2;
  double d2;
  double ad2_x1;
  double ad2_x0;

  double* osBuffer;

  const int npole = 4;
  const int nzero = 4;
  // for samplerate == 48000
  const double* acoeff;
  const double* bcoeff;
  double gain;

  double* xv1;
  double* yv1;
  double* xv2;
  double* yv2;

  enum InputParams { Input, AntiDerivativeLevel, OverSample, NumInputParams };
  enum Outputs { Out1, NumOutputParams };

};

 
}
