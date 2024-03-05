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

    ~Oversampling();

    class OversamplingStage;

    // void init(int factor, int sc_numSamples, float*& upNormalizedTransitionWidths, float*& upStopbandAmplitudesdB,
                                                       // float*& downNormalizedTransitionWidths, float*& downStopbandAmplitudesdB);

    // void addOverSamplingStage(float normalizedTransitionWidthUp, float stopbandAmplitudedBUp,
                              // float normalizedTransitionWidthDown, float stopbandAmplitudedBDown);

    void reset();
    void init(const int& factor, const int& buffer_size);

    void addOverSamplingStage(float normalizedTransitionWidthUp,
                              float stopbandAmplitudeBUp,
                              float normalizedTransitionWidthDown,
                              float stopbandAmplitudeBDown,
                              const int& buffer_size);

    void processSamplesUp(const float* inputBlock, const int& size);
    void processSamplesDown() noexcept;
    std::vector<double> getProcessedSamples();

    // void clearOverSamplingStages();
    

    // int numChannels = 1;
private:

    enum Processed {Up = true, Down = false};
    bool justProcessed{Processed::Up};
    std::vector<OversamplingStage*>* stages;
    int factorOversampling;
    int numStages{0};
};

} // Namespace Oversampling
