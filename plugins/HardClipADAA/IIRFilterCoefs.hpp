// IIRFilterCoefs.hpp
// James Squires (squires.jr@gmail.com)

#pragma once


#include <assert.h>

class BaseIIRFilterCoefs {

public:
  BaseIIRFilterCoefs() = default;
  virtual ~BaseIIRFilterCoefs() {}
  virtual void gen_coefs(int) = 0;
};

class IIRFilterCoefs : public BaseIIRFilterCoefs
  {
  public:
    IIRFilterCoefs() = default;
    virtual ~IIRFilterCoefs() {}

    void gen_coefs(int sampleRate) override {

      assert(sampleRate == SampleRate::FOURTY_FOUR
             || sampleRate == SampleRate::FOURTY_EIGHT
             || sampleRate == SampleRate::EIGHTY_EIGHT);
             // || sampleRate == SampleRate::NINETY_SIX);


      if (sampleRate == SampleRate::FOURTY_FOUR) {
        acoeff = acoeff_41;
        bcoeff = bcoeff_41;
        gain = &gain_41;

      } else if (sampleRate == SampleRate::FOURTY_EIGHT) {
        acoeff = acoeff_48;
        bcoeff = bcoeff_48;
        gain = &gain_48;

      } else if (sampleRate == SampleRate::EIGHTY_EIGHT) {
        acoeff = acoeff_88;
        bcoeff = bcoeff_88;
        gain = &gain_88;

      } 
      /**
       else if (sampleRate == SampleRate::NINETY_SIX) {
        acoeff = acoeff_96;
        bcoeff = bcoeff_96;
        gain = &gain_96;
      }
      **/
    }

    const double* acoeff;
    const double* bcoeff;
    const double* gain;

  private:

    const double acoeff_41[5]={0.017694389856908303,-0.008210396879781378,0.486513017145394,-0.039139137475620284,1};
    const double bcoeff_41[5]={1,4,6,4,1};
    const double gain_41=10.982540095644547;

    // 48k order 4
    const double acoeff_48[5]={0.019451024895140914,0.06459978290184043,0.5150430384020538,0.30290417120437346,1};
    const double bcoeff_48[5]={1,4,6,4,1};
    const double gain_48=8.412206455316431;

    // 88200 order 4
    const double acoeff_88[5]={0.9211819291912352,3.7603495076945235,5.757076379118063,3.9179078653919857,1};
    const double bcoeff_88[5]={1,4,6,4,1};
    const double gain_88=1.0419030157592175;

    // 96000 order 4


    enum SampleRate { FOURTY_FOUR = 44100, FOURTY_EIGHT = 48000, EIGHTY_EIGHT = 88200, NINETY_SIX = 96000 };
    // sample rate == 41000
    // const double acoeff_41[11]={0.000030237345641972724,-0.0004896180481441244,0.005308406339038287,-0.024254250143616325,0.11830637861443252,-0.2604783904269854,0.7506402795579332,-0.8563593626461209,1.6058785470905348,-0.792723694970645,1};
    // const double bcoeff_41[11]={1,10,45,120,210,252,210,120,45,10,1};
    // const double gain_41=662.4150776613988;
    // for samplerate == 48000
    // const double acoeff_48[11] = {0.00002949793572559528,0.00047278455336601714,0.005192457448778824,0.02345936072753368,0.11611377967648272,0.2524443697901383,0.7399674511217426,0.8319489408254346,1.5924340858916823,0.7723958369018249,1};
    // const double bcoeff_48[11] = {1,10,45,120,210,252,210,120,45,10,1};
    // const double gain_48 = 191.95950020926549;

    // sample rate == 88200

    // sample rate == 96000



  };

