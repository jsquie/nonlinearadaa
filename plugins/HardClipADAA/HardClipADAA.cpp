// PluginHardClipADAA.cpp
// James Squires (squires.jr@gmail.com)

#include "SC_PlugIn.hpp"
#include "HardClipADAA.hpp"

static InterfaceTable* ft;

namespace HardClipADAA {

HardClipADAA::HardClipADAA() {
    mCalcFunc = make_calc_function<HardClipADAA, &HardClipADAA::next>();
    next(1);
}

void HardClipADAA::next(int nSamples) {

    // Audio rate input
    const float* input = in(0);

    // Control rate parameter: gain.
    const float gain = in0(1);

    // Output buffer
    float* outbuf = out(0);

    // simple gain function
    for (int i = 0; i < nSamples; ++i) {
        outbuf[i] = input[i] * gain;
    }
}

} // namespace HardClipADAA

PluginLoad(HardClipADAAUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<HardClipADAA::HardClipADAA>(ft, "HardClipADAA", false);
}
