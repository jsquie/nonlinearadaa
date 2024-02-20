// PluginHardClipADAA.hpp
// James Squires (squires.jr@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace HardClipADAA {

class HardClipADAA : public SCUnit {
public:
    HardClipADAA();

    // Destructor
    // ~HardClipADAA();

private:
    // Calc function
    void next(int nSamples);

    // Member variables
};

} // namespace HardClipADAA
