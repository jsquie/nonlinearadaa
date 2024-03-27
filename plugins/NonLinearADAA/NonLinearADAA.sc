HardClipADAA : UGen {

	*ar { |input, adLevel=2, overSample=1|
    if(input.rate!='audio'){input = K2A.ar(input)};
    if((adLevel < 1) && (adLevel > 2)){ Error("adlevel must be either 1 or 2").throw };
    if((overSample < 1) && (overSample > 4)){ Error("Oversample must be a value between 0 and 4").throw };
		^this.multiNew(\audio, input, adLevel, overSample);
	}

	checkInputs {
		^this.checkValidInputs;
	}
}

TanhADAA : UGen {
  *ar { |input, adLevel=2, overSample=1|
    if(input.rate!='audio') {input = K2A.ar(input)};
    ^this.multiNew(\audio, input, adLevel, overSample);
  }

  checkInputs {
    ^this.checkValidInputs;
  }
}

