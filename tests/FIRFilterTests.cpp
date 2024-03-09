#include <gtest/gtest.h>

#include "../plugins/HardClipADAA/FIRFilter.hpp"

TEST(ConvolveTests, BasicConvolve) {
  // convolve takes input sample, circular buffer
  // kernel, kernel size, and buffPosition arguments
  float input[] = {1.0, 2.0};
  float kernel[] = {1.0, 2.0, 3.0};
  float output[] = {0.0, 0.0, 0.0, 0.0};

  ASSERT_EQ(FIRFilter::filterConvolve<float>(input, kernel, output,
                                             unsigned int buffer_size,
                                             unsigned int kernel_size), result)
}
