HardClipADAA : UGen {

	*ar { |input, adLevel=2|
    if(input.rate!='audio'){input = K2A.ar(input)};
		^this.multiNew(\audio, input, adLevel);
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
