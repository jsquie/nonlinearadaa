// PluginHardClipADAA.hpp
// James Squires (squires.jr@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace HardClipADAA {

class HardClipADAA : public SCUnit {

public:
    HardClipADAA();
    // Destructor
    ~HardClipADAA();

private:
    // Calc function
    void next(int nSamples);
    inline void coefs(float* c_buff, const int fc);
    inline void stage_1_filter(int nSamples);
    inline void stage_2_filter(int nSamples);
    static const int nTaps = 26;
    float* coefs_buffer;
    float* s1_prevs_buf;
    float* s2_prevs_buf;
    float* up_sample_buffer;
    float* us_stage1filtered_buf;
    float* us_process_buf;
    float* us_stage2filtered_buf;
    
    int s1_prevs_ptr;
    int s2_prevs_ptr;

    // Member variables
};

} // namespace HardClipADAA
