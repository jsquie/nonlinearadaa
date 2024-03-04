// Limited implementation of oversampling 
// Adapted from JUCE by James Squires
// https://github.com/juce-framework/JUCE/blob/master/modules/juce_dsp/processors/juce_Oversampling.h

#pragma once

#include <memory>

namespace Oversampling
{

 




class Oversampling
{
public:
    Oversampling() = default;

    virtual ~Oversampling() {}

    class OversamplingStage;

    void init(int factor, int sc_numSamples, float*& upNormalizedTransitionWidths, float*& upStopbandAmplitudesdB,
                                                       float*& downNormalizedTransitionWidths, float*& downStopbandAmplitudesdB);

    void addOverSamplingStage(float normalizedTransitionWidthUp, float stopbandAmplitudedBUp,
                              float normalizedTransitionWidthDown, float stopbandAmplitudedBDown);

    /**

    void processSamplesUp(const float* inputBlock) noexcept;
    void processSamplesDown(float*& outputBlock) noexcept;


    void clearOverSamplingStages();
    **/

    // int numChannels = 1;
private:

    int factorOversampling{2};
    int numSamples;
    int stages_add_idx{0};
    std::shared_ptr<std::shared_ptr<OversamplingStage>[]> stages;
};

} // Namespace Oversampling
