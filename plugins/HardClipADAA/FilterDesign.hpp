
#pragma once
#include <memory>
#include <vector>
namespace FilterDesign
{

typedef std::vector<double> Coefficients;

class Filter
{
public:
    Filter() = default;
    ~Filter() {};
};


class IIRLowpassHalfBandPolyphase : public Filter
{
  public:

    IIRLowpassHalfBandPolyphase() {};
    ~IIRLowpassHalfBandPolyphase();

    void designIIRLowpassHalfBandPolyphaseAllpassMethod(const float& normalisedTransitionWidth, const float& stopbandAmplitudedB);
    const int getDirectPathSize() const;
    const int getDelayedPathSize() const;

    std::vector<double> directPath;
    std::vector<double> delayedPath;
    std::vector<double> alpha;

  private:

    int directPath_size;
    int delayedPath_size;

};

} // namespace FilterDesign





