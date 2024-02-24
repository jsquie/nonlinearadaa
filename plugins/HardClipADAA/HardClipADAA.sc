HardClipADAA : UGen {
	*ar { |input, amp = 1, oversample = 1|
		/* TODO */
    if(input.rate!='audio'){input = K2A.ar(input)};
    if(amp.rate!='audio'){amp = K2A.ar(amp)};
		^this.multiNew(\audio, input, amp, oversample);
	}
	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}
