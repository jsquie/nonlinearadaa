HardClipADAA : UGen {

	*ar { |input, adLevel=2, sampleRate=2|
    if(input.rate!='audio'){input = K2A.ar(input)};
		^this.multiNew(\audio, input, adLevel, sampleRate);
	}

	checkInputs {
		^this.checkValidInputs;
	}
}

TanhADAA : UGen {
  *ar { |input, adLevel=2|
    if(input.rate!='audio') {input = K2A.ar(input)};
    ^this.multiNew(\audio, input, adLevel);
  }

  checkInputs {
    ^this.checkValidInputs;
  }
}
