

#include "FilterDesign.hpp"
#include <cmath>
#include <assert.h>
#include "Utils.hpp"
#include <iostream>
#include <memory>

namespace FilterDesign
{


IIRLowpassHalfBandPolyphase::~IIRLowpassHalfBandPolyphase() { std::cout << "Deleting this filter IIRLowpassHalfBandPolyphase object!" << std::endl; };

const int IIRLowpassHalfBandPolyphase::getDirectPathSize() const { return directPath_size; };

const int IIRLowpassHalfBandPolyphase::getDelayedPathSize() const { return delayedPath_size; };

void IIRLowpassHalfBandPolyphase::designIIRLowpassHalfBandPolyphaseAllpassMethod(const float& normalisedTransitionWidth, const float& stopbandAmplitudedB)
{
  // Assert provided transition width and stopband apmplitude dB are within valid ranges
  assert(normalisedTransitionWidth > 0 && normalisedTransitionWidth <= 0.5);
  assert(stopbandAmplitudedB > -300 && stopbandAmplitudedB < -10);

  // std::cout << "Here at the design" << std::endl;

  // Calculate the angular transition width and stopband attenuation as preliminary params for the filter design
  const double wt = Constants::TWO_PI * normalisedTransitionWidth;
  const double ds = 20 * std::log(std::abs(stopbandAmplitudedB));

  // calculate the elliptic function parameters k, kp (the complementary modulus), and elliptic integral modulus e
  double k = std::pow(std::tan ((Constants::PI - wt) / 4), 2.0);
  double kp = std::sqrt (1.0 - k * k);
  double e = (1 - std::sqrt (kp)) / (1 + std::sqrt (kp)) * 0.5;

  // calculate the selectivity factor q based on the elliptic integral modulus e
  double q = e + 2 * std::pow(e, 5.0) + 15 * std::pow(e, 9.0) + 150 * std::pow(e, 13.0);

  // Determine the filter order n using the stopband attenuation and the selectivity factor q
  double k1 = ds * ds / (1 - ds * ds);
  int n = static_cast<int>(int (std::ceil (std::log(k1 * k1 / 16) / std::log(q))));

  // Adjust the filter order n to ensure it is odd, increasing it if necessary
  if (n % 2 == 0) { ++n; }
  if (n == 1) { n = 3; }

  // Calculate the tranformed selectivity factor q1 and a new parameter k1 for the final filter design 
  double q1 = std::pow(q, (double)n);
  k1 = 4 * std::sqrt(q1);

  // Calculate the number of allpass filters needed in the polyphase structure
  const int N = (n - 1) / 2;
  // ai = coefficients

  // std::cout << "Number of coefs N is: " << N << std::endl;
  // Calculate the coefficients for each allpass filter in the structure
  for (int i = 1; i <= N; ++i) {
    double num = 0.0; // Numerator for the Wi calculation
    double delta = 1.0; // Used to check convergence of the infinite series
    int m = 0; // Iterator for the infinite series
    
    // Calculate the numerator of Wi using an infinite series, stopping when delta is sufficiently small
    while (std::abs(delta) > 1e-100) {
      delta = std::pow(-1, m) * std::pow(q, m * (m + 1)) * std::sin((2 * m + 1) * Constants::PI / (double) n);
      num += delta;
      m++;
    }
    num *= 2 * std::pow(q, 0.25);

    double den = 0.0; // Denominator for the Wi calculation
    delta = 1.0; // Reset delta for denominator calculation
    m = 1; // Reset iterator for denominator series

    // Calculate the denominator of Wi in same manner to numerator 
    while (std::abs(delta) > 1e-100) {
      delta = std::pow(-1, m) * std::pow(q, m * m) * std::cos(m * Constants::PI * i / (double)n);
      den += delta;
      ++m;
    }
    den = 1 + 2 * den;
    
    // Calculate the coefficient api and then the allpass coefficient for the ith filter
    double wi = num / den;
    double api = std::sqrt((1 - wi * wi * k) * (1 - wi * wi / k)) / (1 + wi * wi);
    api = (1 - api) / (1 + api);
    std::cout << "Adding coef to alpha: " << api << std::endl;
    alpha.push_back(api);
  }

  
  std::cout << "directPath_size: " << N / 2 << std::endl;
  directPath_size = N / 2;
  delayedPath_size = (N / 2) + 1;

  // Make sure that directpath is sizeof(double * 6) * ((N * 2) + 1) AKA sizeof(Coefficients)
  // Populate the direct path of the polyphase structure with the coefficients calculated above
  
  for (int i = 0; i < N / 2; ++i) {
    double coef_i = alpha[i << 1];
    std::cout << "Adding coefficient to direct path" << std::endl;
    directPath.push_back(coef_i);
  }

  for (int i = 0; i < N / 2; ++i) {
    double coef_i = alpha[(i << 1) + 1];
    std::cout << "Adding coefficient to delayed path" << std::endl;
    delayedPath.push_back(coef_i);
  }



}

}
