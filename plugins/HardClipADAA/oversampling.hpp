// Limited implementation of oversampling
// Adapted from JUCE by James Squires
// https://github.com/juce-framework/JUCE/blob/master/modules/juce_dsp/processors/juce_Oversampling.h

#pragma once

#include "./Filter.hpp"
#include <memory>


namespace JSCDSP::Oversampling {

class Oversampling {
 public:
  Oversampling() = default;

  ~Oversampling();


  void reset();
  void init(const int& factor, const int& buffer_size);

  void addOverSamplingStage(float normalizedTransitionWidthUp,
                            float stopbandAmplitudeBUp,
                            float normalizedTransitionWidthDown,
                            float stopbandAmplitudeBDown);

  void processSamplesUp();
  void processSamplesDown();
  std::vector<double> getProcessedSamples();
  void createFilters();
  Filter::FilterStructure<double>* getFilterStructures() { return filterStructures; };

  // int numChannels = 1;
 private:
  enum Processed { Up = true, Down = false };
  bool justProcessed{Processed::Up};
  int factorOversampling{2};
  int numStages{0};
  Filter::FilterStructure<double>* filterStructures;
};


} // namespace JSCDSP::Oversampling
