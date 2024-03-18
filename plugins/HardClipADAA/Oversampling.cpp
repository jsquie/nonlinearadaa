// Copyright 2024 James Squires
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#include <TargetConditionals.h>
#endif

#include <iostream>

#include "Oversampling.hpp"
#include "OversamplingStage.hpp"

namespace Oversampling {

void Oversampling::init(const int &initOSFactor, const int &initNSamples,
                        const int &M) {
  factor = initOSFactor;
  nSamples = initNSamples;

  for (int i = 0; i < factor - 1; ++i) {
    up_sample_stages.emplace_back(std::ceil(static_cast<double>(M) / 2));
    down_sample_stages.emplace_back(std::ceil(static_cast<double>(M) / 2));
  }
  up_sample_stages.emplace_back(std::floor(static_cast<double>(M) / 2));
  down_sample_stages.emplace_back(std::floor(static_cast<double>(M) / 2));

  for (auto st : up_sample_stages) {
    st.reset();
  }

  for (auto st : down_sample_stages) {
    st.reset();
  }
}

double Oversampling::convolve(const double &input, OversamplingStage &stage,
                              std::shared_ptr<double[]> const &kernel) {
  stage.y0 = 0.0f;
  stage.y1 = 0.0f;

  stage.data[stage.pos] = input;

  int first_inner_prods = stage.size - stage.pos;

#if TARGET_OS_MAC

  vDSP_dotprD(kernel.get() + stage.pos, 1, stage.data.get(), 1, &stage.y0,
              first_inner_prods);
  vDSP_dotprD(kernel.get(), 1, &stage.data.get()[stage.size - stage.pos], 1,
              &stage.y1, stage.pos);

#endif

  stage.pos = (stage.pos == 0) ? stage.size - 1 : stage.pos - 1;
  return stage.y0 + stage.y1;
}

void Oversampling::processSamplesUp(
    const float *const &input, std::vector<std::shared_ptr<double[]>> &kernels,
    std::vector<double> &osBuffer) {
  assert(input != nullptr);
  assert(!kernels.empty());

  for (int n = 0; n < nSamples; ++n) {
    for (int j = 0; j < factor; ++j) {
      if (n == 0) {
        std::cout << "******************* n: " << n << std::endl;
        std::cout << "Conving stage: " << j << std::endl;
        for (int k = 0; k < up_sample_stages.at(j).size - 1; ++k) {
          std::cout << up_sample_stages.at(j).data.get()[k] << ", ";
        }
        std::cout << up_sample_stages.at(j)
                         .data.get()[up_sample_stages.at(j).size - 1]
                  << std::endl;
        std::cout << "with kernel: " << j << std::endl;
        for (int k = 0; k < up_sample_stages.at(j).size - 1; ++k) {
          std::cout << kernels.at(j).get()[k] << ", ";
        }
        std::cout << kernels.at(j).get()[up_sample_stages.at(j).size - 1]
                  << std::endl;
      }
      osBuffer.push_back(convolve(static_cast<double>(input[n]),
                                  up_sample_stages.at(j), kernels.at(j)));
    }
  }
};

void Oversampling::processSamplesDown(
    float *const &output, std::vector<std::shared_ptr<double[]>> &kernels,
    double *const &osBuffer) {
  assert(output != nullptr);
  assert(osBuffer != nullptr);
  assert(!kernels.empty());

  for (int n = 0; n < nSamples; ++n) {
    output[n] = 0.0;

    for (int j = 0; j < factor; ++j) {
      // std::cout << "convolving osBuffer[n << factor + j]: " << ((n << factor)
      // + j) << std::endl;
      output[n] += convolve(osBuffer[(n << 1) + j], down_sample_stages.at(j),
                            kernels.at(j));
    }
  }
}

}  // namespace Oversampling
