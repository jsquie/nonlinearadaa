// Copyright 2024 James Squires

#include "Oversampling.hpp"

#include <assert.h>

#include <algorithm>
#include <iostream>
#include <memory>

#include "./FilterDesign.hpp"

namespace Oversampling {

class Oversampling::OversamplingStage {
 public:
  explicit OversamplingStage(int newFactor) : factor(newFactor) {}

  virtual ~OversamplingStage() {}
  virtual void reset() {}
  virtual void processSamplesUp(const std::vector<double>& inputBlock) {}
  virtual std::vector<double>* getProcessedSamples() { return nullptr; }

  std::vector<double>* buffer;
  std::vector<double>* v1Up;
  int factor;
};

class Oversampling2TimesPolyphaseIIR : public Oversampling::OversamplingStage {
  using ParentType = typename Oversampling::OversamplingStage;

 public:
  Oversampling2TimesPolyphaseIIR(float normalizedTransitionWidthUp,
                                 float stopbandAmplitudeBUp,
                                 float normalizedTransitionWidthDown,
                                 float stopbandAmplitudeBDown, int bufferSize)
      : ParentType(2), buffer_size(bufferSize) {
    aaFilter = std::unique_ptr<FilterDesign::IIRLowpassHalfBandPolyphase>(
        new FilterDesign::IIRLowpassHalfBandPolyphase());
    aiFilter = std::unique_ptr<FilterDesign::IIRLowpassHalfBandPolyphase>(
        new FilterDesign::IIRLowpassHalfBandPolyphase());
    aiFilter->designIIRLowpassHalfBandPolyphaseAllpassMethod(
        normalizedTransitionWidthUp, stopbandAmplitudeBUp);
    aiFilter->designIIRLowpassHalfBandPolyphaseAllpassMethod(
        normalizedTransitionWidthDown, stopbandAmplitudeBDown);

    buffer = new std::vector<double>(buffer_size, 0.0);
    v1Up = new std::vector<double>(
        aiFilter->getDirectPathSize() + aiFilter->getDelayedPathSize(), 0.0);
    assert(buffer != nullptr);
    assert(v1Up != nullptr);
  }

  ~Oversampling2TimesPolyphaseIIR() {
    std::cout << "Deleting 2time poly" << std::endl;
  }

  void reset() override {
    std::fill(buffer->begin(), buffer->end(), 0.0);
    std::fill(v1Up->begin(), v1Up->end(), 0.0);
  }

  /*
  void processSamplesUp(const std::vector<double>& inputBlock) override {
    assert(buffer != nullptr);
    for (double d : inputBlock) {
      buffer->push_back(d);
    }
  }
  */

  void processSamplesUp(const std::vector<double>& inputBlock) override {
    assert(buffer != nullptr);
    assert(v1Up != nullptr);
    assert(aiFilter != nullptr);

    std::cout << "Processing samples Up from inside 2time poly" << std::endl;
    // coefs == aiFilter.directPath then aiFilter.delayedPath
    // int numStages = aiFilter.getDirectPathSize() +
    // aiFilter.getDelayedPathSize();

    // Processing
    // get pointers to the buffer for writing the output, to the delay line, and
    // to the input buffer, v1Up, samples
    int i = 0;

    // loop through each sample
    for (double s : inputBlock) {
      // direct path cascaded allpass filters
      // start with the current input sample
      auto input = s;
      std::cout << "Processing s: " << s
                << " out of input size: " << inputBlock.size() << std::endl;
      std::cout << "Buffersize: " << buffer_size << std::endl;
      std::cout << "aiFilter direct path size: "
                << aiFilter->getDirectPathSize() << std::endl;

      // process the input through each direct path allpass filter
      for (int n = 0; n < aiFilter->getDirectPathSize(); ++n) {
        assert(aiFilter->directPath.empty());
        std::cout << "directPath[n]: " << aiFilter->directPath[n] << std::endl;
        auto alpha = aiFilter->directPath[n];  // the filter coefficient
        std::cout << "v1Up->at(n): " << v1Up->at(n) << std::endl;
        auto output = alpha * input + v1Up->at(n);  // the filtered output
        v1Up->at(n) =
            input - alpha * output;  // update the delay line with the new value
        input = output;              // Pass the output as the next input
      }

      // Output the processed sample to the even-indexed positions in the buffer
      // (because we're upsampling by 2)
      buffer->at(i << 1) = input;

      // Delayed path cascaded allpass filters
      // Reset input to the original sample for processing through the delayed
      // path
      input = s;

      auto delayedPath_offset = aiFilter->getDirectPathSize();

      // process input through each delayed path allpass filter
      for (int n = 0; n < aiFilter->getDelayedPathSize(); ++n) {
        auto alpha = aiFilter->delayedPath[n];  // filter coefficient
        auto output = alpha * input *
                      v1Up->at(n + delayedPath_offset);  // filtered output
        v1Up->at(n + delayedPath_offset) =
            input - alpha * output;  // update the delay line with the new value
        input = output;              // pas the output as teh next input
      }

      // Output the processed sample to the odd-indexed positions in the output
      // buffer
      buffer->at((i << 1) + 1) = input;
      // *buffer[(i << 1) + 1] = input;
      ++i;
    }
    std::cout << "Done processing input block in processing samples up"
              << std::endl;
  }

  std::vector<double>* getProcessedSamples() override {
    assert(buffer != nullptr);
    return buffer;
  }

  std::unique_ptr<FilterDesign::IIRLowpassHalfBandPolyphase> aiFilter;
  std::unique_ptr<FilterDesign::IIRLowpassHalfBandPolyphase> aaFilter;
  const int buffer_size;
};

Oversampling::~Oversampling() {
  std::cout << "Deleting this Oversampling object" << std::endl;
  for (auto s : *stages) {
    delete s;
  }
}

void Oversampling::reset() {
  for (auto stage : *stages) {
    stage->reset();
  }
}

void Oversampling::init(const int& factor, const int& buffer_size) {
  std::cout << "Initializing the oversampling object!" << std::endl;
  factorOversampling = factor;
  stages = new std::vector<OversamplingStage*>;

  for (int n = 0; n < factor; ++n) {
    auto twUp = 0.12f * (n == 0 ? 0.5f : 1.0f);
    auto twDown = 0.15f * (n == 0 ? 0.5f : 1.0f);

    auto gaindBStartUp = -70.0f;
    auto gaindBStartDown = -60.0f;
    auto gaindBFactorUp = 8.0f;
    auto gaindBFactorDown = 8.0f;

    addOverSamplingStage(
        twUp, gaindBStartUp + gaindBFactorUp * static_cast<float>(n), twDown,
        gaindBStartDown + gaindBFactorDown * static_cast<float>(n),
        buffer_size * (n + 1));
  }
}

void Oversampling::addOverSamplingStage(float normalizedTransitionWidthUp,
                                        float stopbandAmplitudeBUp,
                                        float normalizedTransitionWidthDown,
                                        float stopbandAmplitudeBDown,
                                        const int& buffer_size) {
  assert(stages && stages != nullptr);
  std::cout << "Adding a over sampling stage!" << std::endl;
  numStages++;
}

void Oversampling::processSamplesUp(const float* inputBlock, const int& size) {
  assert(stages && stages != nullptr);
  std::cout << "Processing samples Up! From the Oversampling class itself"
            << std::endl;
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
}

void Oversampling::processSamplesDown() noexcept {
  justProcessed = Processed::Down;
}

std::vector<double> Oversampling::getProcessedSamples() {
  // std::cout << "Getting OS processed samples" << std::endl;
  // std::cout << "stages front buffer[0]: " <<
  // (*stages->back()).buffer->front() << std::endl;
  // assert(!(*stages->back()).buffer->empty());
  (justProcessed == Processed::Up)
      ? assert((*stages->back()).buffer != nullptr)
      : assert((*stages->front()).buffer != nullptr);
  return (justProcessed == Processed::Up) ? *(*stages->back()).buffer
                                          : *(*stages->front()).buffer;
}
/**
  void processSamplesDown(double* output, int numSamples)
  {

    assert(bufferSetExternally);
    // loop through each sample
    for (int i = 0; i < numSamples; ++i) {

      // direct path cascaded allpass filters
      auto input = buffer[i << 1];

      // process the input through each direct path allpass filter
      for (int n = 0; n < aaFilter.getDirectPathSize(); ++n) {
        auto alpha = aaFilter.directPath[n].coefficients[0]; // the filter
coefficient auto output = alpha * input + v1Up[n]; // the filtered output
        v1Up[n] = input - alpha * output; // update the delay line with the new
value input = output;  // Pass the output as the next input
      }

      // Output the processed sample to the even-indexed positions in the buffer
(because we're upsampling by 2) auto directOut = input;

      // Delayed path cascaded allpass filters
      input = buffer[(i << 1) + 1];

      auto delayedPath_offset = aaFilter.getDirectPathSize();

      // process input through each delayed path allpass filter
      for (int n = 0; n < aaFilter.getDelayedPathSize(); ++n) {
        auto alpha = aaFilter.delayedPath[n].coefficients[0]; // filter
coefficient auto output = alpha * input * v1Up[n + delayedPath_offset]; //
filtered output v1Up[n + delayedPath_offset] = input - alpha * output; // update
the delay line with the new value input = output; // pas the output as the next
input
      }

      // Output
      output[i] = (delay + directOut) * 0.5;
      delay = input;

    }
  }


**/

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

}  // namespace Oversampling
