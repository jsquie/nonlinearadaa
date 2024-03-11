// Copyright 2024 James Squires

#include "Oversampling.hpp"

#include <assert.h>

#include <iostream>
#include <stdexcept>

namespace JSCDSP::Oversampling {

Oversampling::~Oversampling() {
  std::cout << "Deleting this Oversampling object" << std::endl;
}

void Oversampling::reset() {}

void Oversampling::init(const int& factor, const int& buffer_size) {
  std::cout << "Initializing the oversampling object!" << std::endl;
  factorOversampling = factor;

  for (int n = 0; n < factor; ++n) {
    auto twUp = 0.12f * (n == 0 ? 0.5f : 1.0f);
    auto twDown = 0.15f * (n == 0 ? 0.5f : 1.0f);

    auto gaindBStartUp = -70.0f;
    auto gaindBStartDown = -60.0f;
    auto gaindBFactorUp = 8.0f;
    auto gaindBFactorDown = 8.0f;
    
    // design these filters
    // start with these specs

    /**
    addOverSamplingStage(
        twUp, gaindBStartUp + gaindBFactorUp * static_cast<float>(n), twDown,
        gaindBStartDown + gaindBFactorDown * static_cast<float>(n));
        **/
  }
}

void Oversampling::addOverSamplingStage(float normalizedTransitionWidthUp,
                                        float stopbandAmplitudeBUp,
                                        float normalizedTransitionWidthDown,
                                        float stopbandAmplitudeBDown) {
  std::cout << "Adding a over sampling stage!" << std::endl;
  throw std::logic_error("addOverSamplingStage not yet implemented yet!");
}

void Oversampling::processSamplesUp() {
  std::cout << "Processing samples Up! From the Oversampling class itself"
            << std::endl;
  throw std::logic_error("processSamplesUp not yet implemented yet!");
  /**
  // initial stage uses inputBlock
  // essentially a dumby os stage that does nothing
  for (int i = 0; i < size; ++i) {
    assert((*stages)[0]->buffer != nullptr);
    assert((*stages)[0]->buffer);
    // std::cout << "Adding inputBlock[i]: " << inputBlock[i] << " to
    // stages[0].buffer" << std::endl;

    (*stages)[0]->buffer->push_back(static_cast<double>(inputBlock[i]));
  }

  std::cout
      << "Added inputBlock data to initial stage. On to processing stages up"
      << std::endl;

  for (int i = 1; i < numStages; ++i) {
    assert(stages != nullptr);
    assert(stages);
    auto prev = (*stages)[i - 1]->buffer;
    assert(prev != nullptr);
    assert(prev);
    (*stages)[i]->processSamplesUp(*prev);
  }

  justProcessed = Processed::Up;
       **/
}

void Oversampling::processSamplesDown() {
  justProcessed = Processed::Down;
  throw std::logic_error("processSamplesDown not yet implemented yet!");
}

std::vector<double> Oversampling::getProcessedSamples() {
  // std::cout << "Getting OS processed samples" << std::endl;
  // std::cout << "stages front buffer[0]: " <<
  // (*stages->back()).buffer->front() << std::endl;
  // assert(!(*stages->back()).buffer->empty());
  throw std::logic_error("getProcessedSamples not yet implemented!");
  /**  This does not work for more than 2x oversampling
  (justProcessed == Processed::Up)
      ? assert((*stages->back()).buffer != nullptr)
      : assert((*stages->front()).buffer != nullptr);
  return (justProcessed == Processed::Up) ? *(*stages->back()).buffer
                                          : *(*stages->front()).buffer;
  **/
}


}  // namespace JSCDSP::Oversampling
