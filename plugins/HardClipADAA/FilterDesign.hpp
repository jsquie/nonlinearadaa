
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
    virtual ~Filter() {};
};


class IIRLowpassHalfBandPolyphase : public Filter
{
  public:

    IIRLowpassHalfBandPolyphase() {};
    virtual ~IIRLowpassHalfBandPolyphase() = default;

    void designIIRLowpassHalfBandPolyphaseAllpassMethod(const float& normalisedTransitionWidth, const float& stopbandAmplitudedB);
    const int getDirectPathSize() const;
    const int getDelayedPathSize() const;

    std::shared_ptr<Coefficients[]> directPath;
    std::shared_ptr<Coefficients[]> delayedPath;
    std::shared_ptr<double[]> alpha;

  private:

    int directPath_size;
    int delayedPath_size;

};

} // namespace FilterDesign





