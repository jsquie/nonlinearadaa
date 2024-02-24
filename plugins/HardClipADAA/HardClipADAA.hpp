// PluginHardClipADAA.hpp
// James Squires (squires.jr@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"
#include "VariableOversampling.hpp"
#include <array>

namespace HardClipADAA {

class HardClipADAA : public SCUnit {

public:
    HardClipADAA();
    // Destructor
    ~HardClipADAA();
    VariableOversampling<> oversample;
    float m_oversampling_ratio;

private:
    // Calc function
    float next_os(float sig, float gain);
    void next_aa(int nSamples);
    float aD1_hard_clip(float sig);
    float aD2_hard_clip(float sig);
    float next_adaa(float sig);
    float calcD(float x0, float x1);
    float fallback(float x0, float x2);
    float x1;
    float x2;
    const float TOL = 0.0001;

    enum InputParams { Sig, Amp, OverSample };
    enum Outputs { Out1, NumOutputParams };

    float *osBuffer;

    int m_oversamplingIndex{0};
    
};

} // namespace HardClipADAA
