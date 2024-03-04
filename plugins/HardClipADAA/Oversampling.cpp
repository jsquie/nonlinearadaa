
#include "Oversampling.hpp"
#include "FilterDesign.hpp"
#include "Utils.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>

namespace Oversampling
{


class Oversampling::OversamplingStage
{
public:
  OversamplingStage(int newFactor) : factor (newFactor) {}

  virtual ~OversamplingStage() {}

private:
  double* buffer;
  int factor;

};


class Oversampling2TimesPolyphaseIIR : public Oversampling::OversamplingStage
{

using ParentType = typename Oversampling::OversamplingStage;

public:

  Oversampling2TimesPolyphaseIIR(float normalizedTransitionWidthUp,
                               float stopbandAmplitudeBUp,
                               float normalizedTransitionWidthDown,
                               float stopbandAmplitudeBDown) : ParentType(2), aiFilter(new FilterDesign::IIRLowpassHalfBandPolyphase), aaFilter(new FilterDesign::IIRLowpassHalfBandPolyphase)
  {
    aaFilter->designIIRLowpassHalfBandPolyphaseAllpassMethod (normalizedTransitionWidthUp, stopbandAmplitudeBUp);
    aiFilter->designIIRLowpassHalfBandPolyphaseAllpassMethod (normalizedTransitionWidthDown, stopbandAmplitudeBDown);
    

    // give size info of buffer (nsamples) and set v1Up + v1Down sizes
    v1Up_buffer_size = aiFilter->getDirectPathSize() + aiFilter->getDelayedPathSize();
    v1Down_buffer_size = aaFilter->getDirectPathSize() + aaFilter->getDelayedPathSize();

    bufferSetExternally = false;
    // latency? Would require implementing getcoeffs and getPhaseForFrequency
    // Initializing v1Up ? -> maybe need 

  }

private:
    std::shared_ptr<FilterDesign::IIRLowpassHalfBandPolyphase> aiFilter;
    std::shared_ptr<FilterDesign::IIRLowpassHalfBandPolyphase> aaFilter;
    int v1Up_buffer_size;
    int v1Down_buffer_size;
    bool bufferSetExternally;

};

/**

  void setOSBuffer(double* buffer_ptr, int size) override {
    buffer = buffer_ptr;
    buffer_size = size;
    bufferSetExternally = true;
  }

  void reset() override
  {
    ParentType::reset();
    std::fill(v1Up[0], v1Up[aiFilter.getDirectPathSize() + aiFilter.getDelayedPathSize()], 0.0);
    std::fill(v1Down[0], v1Down[aaFilter.getDirectPathSize() + aaFilter.getDelayedPathSize()], 0.0);
  }

  void processSamplesUp(const double* samples, int numSamples) 
  {

    assert(bufferSetExternally);
    // coefs == aiFilter.directPath then aiFilter.delayedPath 
    // int numStages = aiFilter.getDirectPathSize() + aiFilter.getDelayedPathSize();

    // Processing
    // get pointers to the buffer for writing the output, to the delay line, and to the input
    // buffer, v1Up, samples 

    // loop through each sample
    for (int i = 0; i < numSamples; ++i) {

      // direct path cascaded allpass filters
      // start with the current input sample
      auto input = samples[i];

      // process the input through each direct path allpass filter
      for (int n = 0; n < aiFilter.getDirectPathSize(); ++n) {
        auto alpha = aiFilter.directPath[n].coefficients[0]; // the filter coefficient
        auto output = alpha * input + v1Up[n]; // the filtered output
        v1Up[n] = input - alpha * output; // update the delay line with the new value
        input = output;  // Pass the output as the next input
      }

      // Output the processed sample to the even-indexed positions in the buffer (because we're upsampling by 2)
      buffer[i << 1] = input;

      // Delayed path cascaded allpass filters
      // Reset input to the original sample for processing through the delayed path
      input = samples[i];

      auto delayedPath_offset = aiFilter.getDirectPathSize();

      // process input through each delayed path allpass filter
      for (int n = 0; n < aiFilter.getDelayedPathSize(); ++n) {
        auto alpha = aiFilter.delayedPath[n].coefficients[0]; // filter coefficient
        auto output = alpha * input * v1Up[n + delayedPath_offset]; // filtered output 
        v1Up[n + delayedPath_offset] = input - alpha * output; // update the delay line with the new value
        input = output; // pas the output as teh next input
      }

      // Output the processed sample to the odd-indexed positions in the output buffer
      buffer[(i << 1) + 1] = input;

    }
  }

  void processSamplesDown(double* output, int numSamples)
  {

    assert(bufferSetExternally);
    // loop through each sample
    for (int i = 0; i < numSamples; ++i) {

      // direct path cascaded allpass filters
      auto input = buffer[i << 1];

      // process the input through each direct path allpass filter
      for (int n = 0; n < aaFilter.getDirectPathSize(); ++n) {
        auto alpha = aaFilter.directPath[n].coefficients[0]; // the filter coefficient
        auto output = alpha * input + v1Up[n]; // the filtered output
        v1Up[n] = input - alpha * output; // update the delay line with the new value
        input = output;  // Pass the output as the next input
      }

      // Output the processed sample to the even-indexed positions in the buffer (because we're upsampling by 2)
      auto directOut = input;

      // Delayed path cascaded allpass filters
      input = buffer[(i << 1) + 1];

      auto delayedPath_offset = aaFilter.getDirectPathSize();

      // process input through each delayed path allpass filter
      for (int n = 0; n < aaFilter.getDelayedPathSize(); ++n) {
        auto alpha = aaFilter.delayedPath[n].coefficients[0]; // filter coefficient
        auto output = alpha * input * v1Up[n + delayedPath_offset]; // filtered output 
        v1Up[n + delayedPath_offset] = input - alpha * output; // update the delay line with the new value
        input = output; // pas the output as the next input
      }

      // Output 
      output[i] = (delay + directOut) * 0.5;
      delay = input;

    }
  }

private:

  double* v1Up;
  double* v1Down;
  double delay;

  int v1Up_buffer_size;
  int v1Down_buffer_size;

  FilterDesign::IIRLowpassHalfBandPolyphase aiFilter;
  FilterDesign::IIRLowpassHalfBandPolyphase aaFilter;

  bool bufferSetExternally;


}
**/


void Oversampling::init(int factor, int sc_numSamples, float*& upNormalizedTransitionWidths, float*& upStopbandAmplitudesdB,
                                                       float*& downNormalizedTransitionWidths, float*& downStopbandAmplitudesdB)
{

  // stages_add_idx = 0;

  factorOversampling = factor;
  numSamples = sc_numSamples;

  for (int n = 0; n < factor; ++n) {
    auto twUp = 0.12f * (n == 0 ? 0.5f : 1.0f);
    auto twDown = 0.15f * (n == 0? 0.5f : 1.0f);

    auto gaindBStartUp = -70.0f;
    auto gaindBStartDown = -60.0f;
    auto gaindBFactorUp = 8.0f;
    auto gaindBFactorDown = 8.0f;

    upNormalizedTransitionWidths[n] = twUp;
    upStopbandAmplitudesdB[n] = gaindBStartUp + gaindBFactorUp * (float)n;
    downNormalizedTransitionWidths[n] = twDown;
    downStopbandAmplitudesdB[n] = gaindBStartDown + gaindBFactorDown * (float)n;

    // stages[n] = Oversampling2TimesPolyphaseIIR(1.0, 1.0, 1.0, 1.0);
  }

};


void Oversampling::addOverSamplingStage(float normalizedTransitionWidthUp, float stopbandAmplitudedBUp,
                                        float normalizedTransitionWidthDown, float stopbandAmplitudedBDown)
{

  stages[stages_add_idx] = std::shared_ptr<Oversampling2TimesPolyphaseIIR>(new Oversampling2TimesPolyphaseIIR(normalizedTransitionWidthUp, stopbandAmplitudedBUp,
                                                              normalizedTransitionWidthDown, stopbandAmplitudedBDown));
  
  stages_add_idx += 1;
  factorOversampling *= 2;
};

/**
void Oversampling::processSamplesUp(double* samples)
{
  assert(stages != NULL);

  // up sample stage 1 (input)
  auto current_stage = stages[0];

  current_stage->processSamplesUp(samples);

  for (int i = 1; i < stages_add_idx - 1; ++i) {
    auto prev_stage_data = current_stage.buffer;
    current_stage = stages[i];
    current_stage->processSamplesUp(prev_stage_data);
  }

}

void Oversampling::processSamplesDown(double* output)
{
  assert(stages != NULL);

  for (int n = stages_add_idx - 1; n > 0; --n) {
    auto& stage = *stages[n];
    stage.processSamplesDown(stages[n-1].buffer, currNumSamples);
  }

  stages[0].processSamplesDown(output, currNumSamples);

}
**/
} // Namespace Oversampling 
