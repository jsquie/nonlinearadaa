// PluginHardClipADAA.cpp
// James Squires (squires.jr@gmail.com)

#include "SC_PlugIn.hpp"
#include "HardClipADAA.hpp"
#include "SC_PlugIn.h"

static InterfaceTable* ft;

namespace HardClipADAA {

HardClipADAA::HardClipADAA() {
  const float samplerate = (float) sampleRate();

  oversample.reset(samplerate);
  m_oversamplingIndex = sc_clip((int)in0(OverSample), 0, 4);
  oversample.setOversamplingIndex(m_oversamplingIndex);
  osBuffer = oversample.getOSBuffer();

  x1 = 0.0;
  x2 = 0.0;

  mCalcFunc = make_calc_function<HardClipADAA, &HardClipADAA::next_aa>();
  next_aa(1);
}

HardClipADAA::~HardClipADAA() {}

float HardClipADAA::aD1_hard_clip(float sig) {
  if (std::abs(sig) >= 1.0) {
    return 0.5 * (sig * sig);
  } else {
    if (sig < 0) {
      return (sig * -1.0) - 0.5;
    } else {
      return sig - 0.5;
    }
  }
}

float HardClipADAA::aD2_hard_clip(float sig) {
  if (std::abs(sig) <= 1.0) {
    return (1.0 / 6.0) * (sig * sig * sig);
  } else {
    float sign = sig <= 0.0 ? -1.0 : 1.0;
    return (0.5 * (sig * sig) + (1.0 / 6.0)) * sign - (sig * 0.5);
  }
}

float HardClipADAA::next_adaa(float sig) {
  float res;

  if (std::abs(sig - x1) < TOL) {
    res = fallback(sig, x2); 
  } else {
    res = (2.0 / (sig - x2)) * (calcD(sig, x1) - calcD(x1, x2)); 
  }
  x2 = x1;
  x1 = sig;

  return res;
}

float HardClipADAA::calcD(float x0, float x1) {
  if (std::abs(x0 - x1) < TOL) {
    return aD1_hard_clip((x0 + x1) * 0.5);
  } else {
    return (aD2_hard_clip(x0) - aD2_hard_clip(x1))/ (x0 - x1);
  }
}

float HardClipADAA::fallback(float x0, float x2) {
  float x_bar = (x0 + x2) * 0.5;
  float delta = x_bar - x0;

  if (delta < TOL) {
    return sc_clip((x_bar + x0) * 0.5, -1.0, 1.0);
  } else {
    return (2.0 / delta) * (aD1_hard_clip(x_bar) + (aD2_hard_clip(x0) - aD2_hard_clip(x_bar)) / delta);
  }
}

float HardClipADAA::next_os(float sig, float amp) {
  float out;

  oversample.upsample(sig);

  for (int k = 0; k < m_oversampling_ratio; k++) {
    osBuffer[k] = next_adaa(osBuffer[k]);
  }
  if (m_oversamplingIndex != 0) {
    out = oversample.downsample();
  } else {
    out = osBuffer[0];
  }
  return out;
}

void HardClipADAA::next_aa(int nSamples) {
  const float *sig = in(Sig);
  const float *amp = in(Amp);

  float *outbuf = out(Out1);

  for (int i = 0; i < nSamples; ++i) {
    outbuf[i] = next_os(sig[i], amp[i]);
  }
}

} // namespace HardClipADAA

PluginLoad(HardClipADAAUGens) {
  // Plugin magic
  ft = inTable;
  registerUnit<HardClipADAA::HardClipADAA>(ft, "HardClipADAA", true);
}
