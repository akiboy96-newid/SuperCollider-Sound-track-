//Final
//Aakash Kedia

//Track name: Mutated

(
var score, synthDef;
var additiveFunc;
var partialGains = [1.0, 1.67, 1.5, 1.8, 2.67, 1.67, 1.46, 1.33, 1.33, 1.0, 1.33].ampdb; // gains for Risset Bell
var partialRatios = [0.56, 0.56, 0.92, 0.92, 1.19, 1.7, 2.0, 2.74, 3.0, 3.76, 4.07]; // ratios for each partial
var partialDurs = [1.0, 0.9, 0.65, 0.55, 0.325, 0.35, 0.25, 0.2, 0.15, 0.1, 0.075]; // duration scaling;
var partialData, noteArray, rhythmArray;
var sec1Func, sec2Func, sec3Func;
var synthDefBells;

var synthDefShep, additiveFuncShep, partialSizeShep, partialGainsShep = [0];
var partialDataShep, partialRatiosShep = [1.0]; // tuning ratios for each partial
var noteParamsFunc;
var cmRatioFunc, removeDuplicatesFunc, fmRatiosFunc;
var noteNRatios, notePNumbers, noteModIndexs;
var numNotes;
var noteParams, noteParams1,noteParamsMid, start, dur, gain, freq;
var synthDefFM;


//Gran Synth Var
var sgsSynthDef, rmblnSynth, synthDefGran;
var noteParamsGran, startGran, durGran, risGran, decGran, gainGran;
var noteFreqsGran, noteFormFreqsGran, noteQsGran;
var noteParamsFuncGran, additiveFuncGran;
var numOctDivisionsGran, numNotesGran;
var formDataGran, formData1Gran, formData2Gran, formSizesGran;


// fundamental is the starting frequency
// numOctaves is how many octaves we want to glissando through
// startDelay is how many seconds before the next glissando starts and also the time it takes to glissando one octave
// centerFreq is where in the glissando is most of the energy
var numCycles = 25, fundamental = 20, numOctaves = 5, startDelay =9.81, centerFreq = 110, octave = 2;
var	outputPath, headerFormat, sampleFormat, numOutputChannels;


// helper functions
var sinCosPanLaw;



// set the NRT vars here...
outputPath = "~/Desktop/Kedia_Final.wav"; // output file path
headerFormat = "WAV";                 // soundfile header format
sampleFormat = "int24";               // soundfile sample format
numOutputChannels = 2;                // stereo --> 2 channels

// create a score
score = CtkScore.new;

// sine-cosine panning law coefficient function
// angle argument in degrees
sinCosPanLaw = { arg angleInDegrees = 0;
    var angleInRadians;
    var theta;

    angleInRadians = angleInDegrees/180*pi;

    theta = pi/4 - angleInRadians;

    [theta.cos, theta.sin]
};


//GranSynthCode

sgsSynthDef = CtkSynthDef.new(\myGrainSinSynth, {arg dur, gain, ris = 0.1, dec = 0.1, freq = 440.0, formFreq = 1760.0, q = 1.0,  panAngle = 25.0;

	// variables
	var bus = 0;      // var to specify output bus: first output
	var trigger;
	var granSig;
	var out;          // output!
	var amp;          // a few vars for synthesis
	var grainDur, grainFreq, envFreq, wavFreq;
	var ampEnv;       // var for envelope signal


	// initial calcs
	amp = gain.dbamp; // convert from gain in dB to linear amplitude scale
	grainFreq = freq;
	envFreq = q.reciprocal * (formFreq/2);
	wavFreq = formFreq;
	grainDur = envFreq.reciprocal;

	// the amplitude envelope nested in the UGen that synthesises the envelope
	ampEnv = EnvGen.kr(
		Env.linen(ris, 1.0 - (ris + dec), dec),
		timeScale: dur
	);

	// granular (grain frequency) trigger
	trigger = Impulse.ar(grainFreq);

	// granular synthesis
	granSig = GrainSin.ar(trigger: trigger, dur: grainDur, freq: wavFreq);

	// apply the amplitude envelope
	granSig = amp * ampEnv * granSig;

	// expand to two channels - panning
	out = sinCosPanLaw.value(panAngle) * granSig;  // <-- Panning happens here!

	// out!!
	Out.ar(bus, out)
});


// Ring Rodulated (RM) Band-limited Noise (BLN) Synthesis: Noise Band Synthesis
// NOTE: this synthDef is named rmblnSynth
rmblnSynth = CtkSynthDef.new(\myrmblnSynth, {arg dur, gain, ris = 0.1, dec = 0.1, formFreq = 1760.0, q = 1.0,  panAngle = 10.0;

	// variables
	var bus = 0;      // var to specify output bus: first output
	var noise, carOsc;
	var out;          // output!
	var amp;          // a few vars for synthesis
	var carFreq, noiFreq;
	var ampEnv;       // var for envelope signal


	// initial calcs
	amp = gain.dbamp;         // convert from gain in dB to linear amplitude scale
	carFreq = formFreq;       // assign carrier frequency to formFreq
	noiFreq = carFreq/(2*q);  // calculate noiFreq

	// the amplitude envelope nested in the UGen that synthesises the envelope
	ampEnv = EnvGen.kr(
		Env.linen(ris, 1.0 - (ris + dec), dec),
		timeScale: dur
	);

	noise = LFNoise2.ar(noiFreq);              // (quadratic interpolation sample-and-hold noise)
	carOsc = SinOsc.ar(carFreq);               // simple carrier (single sinusoid)

	// apply the amplitude envelope and Ring Modulate
	noise = amp * ampEnv * noise * carOsc;

	// expand to two channels - panning
	out = sinCosPanLaw.value(panAngle) * noise;  // <-- Panning happens here!

	// out!!
	Out.ar(bus, out)
});





synthDefShep = CtkSynthDef.new(\mySinOscSynth, {arg dur, gain = -12, cFreq = 442, startFreq = 110.0, endFreq = 440, panAngle = 5.0;
	var envShep, envGenShep, envTimeShep, ampShep;
	ampShep = gain.dbamp;
	envTimeShep = startDelay * (cFreq/fundamental).log2; // since we have startDelay seconds between octaves, we can multiply by how many octaves is takes to get to the center frequency
	envShep = Env(
		[0, 1, 0], // amp values
		[envTimeShep,  // time to centerFreq
			dur - envTimeShep], // the remainder of the time after centerFreq
		[5, -5] // curves
	);
	envGenShep = EnvGen.kr(envShep);
	ampShep = envGenShep * ampShep;
	Out.ar(0,
		// 'amp' is a value between 0 and 1.
		SinOsc.ar(XLine.ar(startFreq, endFreq, dur), 0, ampShep) * sinCosPanLaw.value(panAngle); // XLine for exponential between freqs

	)
});

// function to add individual notes to our score for additive synthesis
additiveFuncShep = {arg startShep, durShep, gainShep = -14, startFreq = 142.0, endFreq = 442.0, partialDataShep; // pass in values

	// iterate through the partialData array to generate individual notes to add to the score

	partialDataShep.do({arg thisPartialDataShep, iShep;
		var thisPartialFreqShep;
		var thisPartialGainShep;
		var thisPartialRatioShep;

		// retreive partial dur, gain and ratio
		thisPartialGainShep = thisPartialDataShep.at(0);
		thisPartialRatioShep = thisPartialDataShep.at(1);

		thisPartialGainShep= gainShep + thisPartialGainShep; // add gains

		// create a note for each partial...
		score.add(
			synthDefShep.note(startShep, durShep)
			.dur_(durShep)
			.gain_(thisPartialGainShep)
			.startFreq_(startFreq * thisPartialRatioShep)
			.endFreq_(endFreq * thisPartialRatioShep)
			.cFreq_(centerFreq);
		);
	})
};



// group partial durs, gains and ratios into a new array
partialDataShep = [partialGainsShep, partialRatiosShep];
partialDataShep = partialDataShep.flop;


// evaluate the additive synthesis function
// args: start, dur, amp, ris, dec, freq, partialData
// this function adds individual partials to the score to be played

numCycles.do{arg iShep;

	additiveFuncShep.value(
		45+startDelay + iShep, // start times
		startDelay * numOctaves + (iShep), // duration
		-18, // amp
		fundamental * iShep, //startFreq
		fundamental * (octave**numOctaves), // endFreq
		partialDataShep); // harmonic data
};



synthDefBells = CtkSynthDef.new(\mySinOscSynth2, {arg dur, gain = -8, freq = 440.0, ris = 0.1, pan;
	var env, envGen, amp;
	env = EnvGen.kr(Env([0, 1, 0], [ris, dur - ris], [4, -4])); // env define within synthDef
	amp = env * gain.dbamp;
	Out.ar(0, Pan2.ar(
		// 'amp' is a value between 0 and 1.
		SinOsc.ar(freq, 0, amp), pan
	)
	)
});

additiveFunc = {arg start, dur, gain = -12, ris, freq, partialData, pan; // pass in values

	partialData.do({arg thisPartialData, i;
		var thisPartialDur;
		var thisPartialFreq;
		var thisPartialGain;
		var thisPartialRatio;

		// retreive partial dur, gain and ratio
		thisPartialDur = thisPartialData.at(0);
		thisPartialGain = thisPartialData.at(1);
		thisPartialRatio = thisPartialData.at(2);

		thisPartialDur = dur * thisPartialDur; // scale partial dur by dur argument
		thisPartialGain = gain + thisPartialGain; // add gains together
		thisPartialFreq = freq * thisPartialRatio; // multiply freq by index (harmonic series!)


		if (i == 1, //test value 1, if it is second harmonic
			{ thisPartialFreq = thisPartialFreq + 1 }, // function 1

			{if ( i == 3, //test value 2, if it is fourth harmonic
				{ thisPartialFreq = thisPartialFreq + 1.5} // function 2
			)
			}
		);

		// then and add note for each partial to the score
		score.add(
			synthDefBells.note(start+0.1.rrand(1.0), thisPartialDur)
			.dur_(thisPartialDur)
			.gain_(thisPartialGain)
			.freq_(thisPartialFreq)
			.ris_(ris)
			.pan_(pan)
		);

	})
};

 partialData = [partialDurs, partialGains, partialRatios].flop;

rhythmArray = [0.37, 0.47, 0.764, 0.96, 1.25];
// create my note array with rhythms, a random duration and a random amp
noteArray = rhythmArray.size.collect{arg i;
	[rhythmArray[i], 1.5.rrand(2.5), -24.rrand(-18)];
};


// function for the first section of note calls
sec1Func = {arg start, noteArray, fundamental, partialData; // arguments are start time and a note array [start times, durs, amps]
	noteArray.do{arg note, i; //iterate over that array
		var numHarms, harmArray, partialArray, noteDur, noteStart, noteGain;
		noteStart = note[0]; // unpack note start time
		noteDur = note[1]; // unpack note duration
		noteGain = note[2]; // unpack note gain
		numHarms = 3.rrand(5); // pick a random subset of harmonics
		harmArray = Array.fill(numHarms, {0.rrand(partialData.size - 1)}); // fill Array with partial numbers
		// collect into a new array the partial data only at subset harmonics
		partialArray = harmArray.collect{arg index;
			partialData[index];
		};
		additiveFunc.value(
			start + noteStart,  // start time and offset note start time
			noteDur, // note duration
			noteGain,  // note gain
			0.01,  // attack time
			fundamental, // the fundamental frequency of our bell
			partialArray, // our new partialData array with a subset of partials
			-1.0.rrand(1.0);
		);
	};
};

sec1Func.value(0.0, noteArray, 440.0, partialData); // call the function
sec1Func.value(3.0, noteArray, 220.0, partialData); // call the function
sec1Func.value(6.0, noteArray, 110.0, partialData);
sec1Func.value(10.0, noteArray, 330.0, partialData);
sec1Func.value(13.0, noteArray, 330.0, partialData);
sec1Func.value(32.0, noteArray, 220.0, partialData);
sec1Func.value(38.0, noteArray, 330.0, partialData);
sec1Func.value(37.0, noteArray, 170.0, partialData);
sec1Func.value(150.0, noteArray, 440.0, partialData); // call the function
sec1Func.value(183.0, noteArray, 220.0, partialData); // call the function
sec1Func.value(176.0, noteArray, 220.0, partialData);
sec1Func.value(170.0, noteArray, 330.0, partialData);
sec1Func.value(169.0, noteArray, 380.0, partialData);
sec1Func.value(175.0, noteArray, 550.0, partialData);
sec1Func.value(172.0, noteArray, 330.0, partialData);
sec1Func.value(177.0, noteArray, 170.0, partialData);


// second section function to add the slow attack bell
sec2Func = {arg start, dur, fundamental, partialData;
	partialData = partialData.flop;
	partialData = partialData.put(0, Array.fill(partialData[0].size, {1.0})).flop; // put all 1.0s for durations
	additiveFunc.value(
		start, // starttime
		dur, //duration
		-22, //amplitude
		dur/2, // attack time
		fundamental,  // fundamental freq
		partialData,
		0
	);


};


// call sec2Func
sec2Func.value(15, 3.0, 440.0, partialData);
sec2Func.value(36.0, 2.0, 440.0, partialData);
sec2Func.value(30.0, 5.0, 330.0, partialData);
sec2Func.value(25.0, 3.0, 550.0, partialData);
sec2Func.value (23.0, 3.7, 660.0, partialData);
sec2Func.value(127.0, 1.0, 770.0, partialData);
sec2Func.value(125, 3.0, 440.0, partialData);
sec2Func.value(136.0, 2.0, 440.0, partialData);
sec2Func.value(130.0, 5.0, 330.0, partialData);
sec2Func.value(135.0, 3.0, 550.0, partialData);
sec2Func.value (153.0, 3.7, 660.0, partialData);
sec2Func.value(167.0, 1.0, 770.0, partialData);


// add final bell attack
additiveFunc.value(38, 8, -24, 0.01, 440, partialData, 6);
additiveFunc.value(18, 11, -24, 0.05, 880, partialData, -6);
additiveFunc.value(28, 13, -24, 0.05, 330, partialData, 5);
additiveFunc.value(33, 19, -22, 2, 1110, partialData, -5);
additiveFunc.value(39, 20, -26, 0.05, 770, partialData,1);
additiveFunc.value(35, 20, -24, 0.05, 220, partialData, 0);
additiveFunc.value(138, 8, -24, 0.01, 440, partialData, 6);
additiveFunc.value(138, 11, -24, 0.05, 220, partialData, -6);
additiveFunc.value(128, 13, -24, 0.05, 330, partialData, 5);
additiveFunc.value(133, 19, -22, 2, 1110, partialData, -5);
additiveFunc.value(139, 20, -26, 0.05, 770, partialData,1);
additiveFunc.value(135, 20, -24, 0.05, 220, partialData, 0);


synthDef = CtkSynthDef.new(\myCSSynth, {arg dur, gain, ris = 0.1, dec = 0.1, freq = 440.0, numharm = 1000, panAngle = 10.0;

    // variables
    var bus = 0;      // var to specify output bus: first output
    var blip, out;    // vars assigned to audio signals
    var amp, osc;          // a few vars for synthesis
    var ampEnv, env;     // var for envelope signal
	var levels, times, curves; // vars for envelope

	levels = [0, 1, 0.5, 0.25, 0]; // envelope amplitude scales (may be easier to spec in dB)
        times = [0.1, 0.1, 0.9, 0.1]; // envelope times as percentages
        curves = [\lin, \exp, -10, \sin]; // envelope segment curves

    // calcs
    amp = gain.dbamp; // convert from gain in dB to linear amplitude scale

    // the amplitude envelope nested in the UGen that synthesises the envelope
        env = Env.new(levels, times, curves);

        // the UGen that synthesises the envelope
        ampEnv = EnvGen.kr(env, timeScale: dur);

    // complex sources synthesis
    blip = Blip.ar(freq, numharm);

    // envelope
    out = blip * (amp * ampEnv);

    // expand to two channels - panning
    out = sinCosPanLaw.value(panAngle) * out;  // <-- Panning happens here!

    // out!!
    Out.ar(bus, out)
});


// // Example 2
// //
// // Band Limited Impulse -- 1/8 bandwidth
 Array.geom(7, 220.0, 0.5).do({ arg freq, i;
        var dur = 1.7; // duration in seconds
         var gain = -9;
         var nyqFreq = s.sampleRate/2;

         score.add(synthDef.note(starttime: i + dur, duration: dur).dur_(dur).gain_(gain).freq_(freq).numharm_((nyqFreq/8 * freq.reciprocal).asInteger));
 });
 Array.geom(15, 220.0, 1.2).do({ arg freq, i;
        var dur = 2.2; // duration in seconds
         var gain = -9;
         var nyqFreq = s.sampleRate/2;

         score.add(synthDef.note(starttime: i +8, duration: dur).dur_(dur).gain_(gain).freq_(freq).numharm_((nyqFreq/8 * freq.reciprocal).asInteger));
 });

 Array.geom(10, 330.0, 0.8).do({ arg freq, i;
        var dur = 1.2; // duration in seconds
         var gain = -9;
         var nyqFreq = s.sampleRate/2;
score.add(synthDef.note(starttime: i +15, duration: dur).dur_(dur).gain_(gain).freq_(freq).numharm_((nyqFreq/8 * freq.reciprocal).asInteger));
 });

Array.geom(18, 110.0, 0.9).do({ arg freq, i;
        var dur = 0.6; // duration in seconds
         var gain = -9;
         var nyqFreq = s.sampleRate/2;
score.add(synthDef.note(starttime: i +23, duration: dur+i).dur_(dur).gain_(gain).freq_(freq).numharm_((nyqFreq/8 * freq.reciprocal).asInteger));
 });


synthDefFM = CtkSynthDef.new(\myFMSynth, {arg dur, gain, ris = 0.1, dec = 0.1, freq = 440.0, carRatio = 1, modRatio = 1, modIndex = 1.0,  panAngle = 10.0;

    // variables
    var bus = 0;      // var to specify output bus: first output
    var carOsc, modOsc;  // oscillators
    var out;          // output!
    var amp;          // a few vars for synthesis
    var carFreq, modFreq;
    var modDev;
    var ampEnv;       // var for envelope signal


    // initial calcs
    amp = gain.dbamp; // convert from gain in dB to linear amplitude scale
    carFreq = carRatio * freq;
    modFreq = modRatio * freq;
    modDev = modIndex * modFreq;

    // the amplitude envelope nested in the UGen that synthesises the envelope
    ampEnv = EnvGen.kr(
        Env.linen(ris, 1.0 - (ris + dec), dec),
        timeScale: dur
    );

    modOsc = SinOsc.ar(modFreq, 0, modDev);         // simple modulator (single sinusoid)
    carOsc = SinOsc.ar(carFreq + modOsc, 0, amp);   // simple carrier (single sinusoid)

    // apply the amplitude envelope
    carOsc = ampEnv * carOsc;

    // expand to two channels - panning
    out = sinCosPanLaw.value(panAngle) * carOsc;  // <-- Panning happens here!

    // out!!
    Out.ar(bus, out)
});


// note parameter function - function to add notes to score
noteParamsFunc = { arg myParams;

    // construct score - iterate through noteParams array
    myParams.do({arg params;
        score.add(
            synthDefFM.note(
                starttime: params.at(0),                      // starttime
                duration: params.at(1)                        // dur
            ).dur_(params.at(1))                              // dur
            .gain_(params.at(2))                              // gain
            .freq_(params.at(3))                              // freq
            .carRatio_(cmRatioFunc.value(params.at(4), params.at(5)).at(0)) // carRatio
            .modRatio_(cmRatioFunc.value(params.at(4), params.at(5)).at(1)) // modRatio
            .modIndex_(params.at(6))                          // modIndex
			.panAngle_(params.at(7))
        );
    });
};


// function to remove duplicates
removeDuplicatesFunc = { arg array;
    var result;

    result = Array.newClear;

    array.do({ arg item;
        result.includes(item).not.if({
            result = result.add(item);
        })
    });
    result
};


// function to calculate c:m
// n, n >=1 ratio between f0 and f1
// p, p = 0, 1, 2, 3, ... partial number
cmRatioFunc = {arg n, p;
    var cm;

    p.even.if({
        cm = [                        // p is even case
            (p)/2 * (1 + n) + 1,      // carrier
            (1 + n)                   // modulator
        ]
    },{
        cm = [                        // p is odd case
            (p + 1)/2 * (1 + n) - 1,  // carrier
            (1 + n)                   // modulator
            ]
    });
    cm
};


// function to calculate FM spectrum ratios
fmRatiosFunc = {arg cm, k = 4;
    var ratios;

    ratios = (k + 1).collect({arg kNum;
        [cm.at(0) - (kNum * cm.at(1)), cm.at(0) + (kNum * cm.at(1))];
    });

    ratios = ratios.flatten.abs.sort;
    ratios = removeDuplicatesFunc.value(ratios);

    ratios;
};

// note parameter function - function to parse note parameter list...
// ... calls additiveFunc to create score
noteParamsFuncGran = { arg mySynthDefGran, myParams, formDataGran;

	// iterate through noteParams array - call additiveFunc
	myParams.do({arg params;
		additiveFuncGran.value(
			mySynthDefGran,          // mySynthDef
			params.at(0),        // starttime
			params.at(1),        // dur
			params.at(2),        // gain
			params.at(3),        // ris
			params.at(4),        // dec
			params.at(5),        // freq
			formDataGran             // formData
		);
	});
};


// function to add individual notes to our score for additive synthesis
additiveFuncGran = {arg mySynthDefGran, startGran, durGran, gain, risGran = 0.1, decGran = 0.1, freqGran = 440.0, formDataGran; // pass in values

	// iterate through the partialData array to generate individual notes to add to the score
	formDataGran.do({arg thisFormDataGran, iGran;
		var thisFormFreqGran, thisFormGainGran, thisFormQGran;

		// retreive formant freq, gain, Q
		thisFormFreqGran = thisFormDataGran.at(0);
		thisFormGainGran = thisFormDataGran.at(1);
		thisFormQGran = thisFormDataGran.at(2);

		thisFormFreqGran = CtkControl.env(Env(thisFormFreqGran, [1.0, 1.0, 1.0].normalizeSum, \exp), timeScale: durGran);
		thisFormGainGran = CtkControl.env(Env(thisFormGainGran + gain, [1.0, 1.0, 1.0].normalizeSum), timeScale: durGran);
		thisFormQGran = CtkControl.env(Env(thisFormQGran, [1.0, 1.0, 1.0].normalizeSum, \exp), timeScale: durGran);

		// create a note for each formant...
		if ( mySynthDefGran == sgsSynthDef,   // check if tonal or noise
			{
				score.add(                       // <--- for TONAL formants
					mySynthDefGran.note(             // mySynthDef
						starttime: startGran,        // start
						duration: durGran            // dur
					).dur_(durGran)                  // dur
					.gain_(thisFormGainGran)  // gain
					.ris_(risGran)                   // ris
					.dec_(decGran)                   // dec
					.formFreq_(thisFormFreqGran)     // formFreq
					.q_(thisFormQGran)                  // q
				)
			},
			{
				score.add(
					mySynthDefGran.note(             // mySynthDef
						starttime: startGran,        // start
						duration: durGran            // dur
					).dur_(durGran)                  // dur
					.gain_(thisFormGainGran)  // gain
					.ris_(risGran)                   // ris
					.dec_(decGran)                   // dec
					.formFreq_(thisFormFreqGran)     // formFreq
					.q_(thisFormQGran)           // q
				)
			}
		)
	})
};

// -------------------------------------------
// parameters...


// specify fm spectrum -- n = 2, (f1/f0) octave relationship
// n    : f1/f0 ratio
// p    : fc partial number
numNotes = 95;
start = 48;
dur = 0.32;
gain = -24;
freq = [220.0, 440.0, 110.0, 880, 1110, 330, 550, 660, 770, 110,880,440];
noteNRatios = Array.fill(numNotes, 1.65);    // octave between f0 and f1
notePNumbers = Array.fill(numNotes, {0.rrand(4)}); // <--- This is our "data of interest"
noteModIndexs = Array.fill(numNotes, {0.2.rrand(1.64)});

// -------------------------------------------
// use -collect (iteration!) to pack into array, noteParams
noteParams = numNotes.collect({arg i;
    Array.with(
		        start + (1.0.rrand(100.0) * dur),   // start
		        0.5.rrand(dur),                 // dur
		        -36.rrand(gain),                // gain
        freq.choose,                // freq
        noteNRatios.choose.postln,   // nRatio
        notePNumbers.choose,  // pNumber
        noteModIndexs.choose, // index of modulation
		-45.rrand(45)
	)
});

noteParams1 = numNotes.collect({arg i;
    Array.with(
		        100 + start + (1.0.rrand(100.0) * dur),   // start
		        0.5.rrand(dur),                 // dur
		        -36.rrand(gain),                // gain
        freq.choose,                // freq
        noteNRatios.choose.postln,   // nRatio
        notePNumbers.choose,  // pNumber
        noteModIndexs.choose, // index of modulation
		-45.rrand(45)
	)
});

noteParamsMid = numNotes.collect({arg i;
    Array.with(
		        60 + start + (1.0.rrand(100.0) * dur),   // start
		        0.4.rrand(dur),                 // dur
		        -36.rrand(gain),                // gain
        freq.choose,                // freq
        noteNRatios.choose.postln,   // nRatio
        notePNumbers.choose,  // pNumber
        noteModIndexs.choose, // index of modulation
		-45.rrand(45)
	)
});

	// -------------------------------------------
// parameters...

// Formant Data
// [[ formFreq, formGain, formQ ], ... ]

// Example 1
// Male Sung Vowel A - D&J p230
formData1Gran = [
	[  609.0,   0.0, 2.00 ], // 1st formant
	[ 1000.0,  -6.0, 3.00 ], // 2nd formant
	[ 2450.0, -12.0, 5.00 ], // 3rd formant
	[ 2700.0, -11.0, 5.00 ], // 4th formant
	[ 3240.0, -24.0, 5.00 ]  // 5th formant
];

formData2Gran = [
	[  238.0,   0.0, 0.70 ], // 1st formant
	[ 1741.0, -20.0, 3.00 ], // 2nd formant
	[ 2450.0, -16.0, 5.00 ], // 3rd formant
	[ 2900.0, -20.0, 5.00 ], // 4th formant
	[ 4000.0, -32.0, 3.00 ]  // 5th formant
];

formDataGran = [formData1Gran, formData1Gran, formData2Gran, formData2Gran];

formSizesGran = formDataGran.shape;
formDataGran = formDataGran.lace.lace.reshape(formSizesGran.at(2), formSizesGran.at(1), formSizesGran.at(0))
.lace.reshape(formSizesGran.at(1), formSizesGran.at(2), formSizesGran.at(0));

startGran = 140;
durGran = 1.6;
gain = -14;
risGran = 0.1;
decGran = 0.1;
numOctDivisionsGran = 5;
numNotesGran = numOctDivisionsGran + 1; // ascend a whole octave
noteFreqsGran = 440.0 * Array.geom(numNotesGran, 1, 2.pow(1/(numOctDivisionsGran)));
synthDefGran = sgsSynthDef;   // <--- This is our "data of interest", choose TONAL formants
// synthDef = rmblnSynth;      // <--- This is our "data of interest", choose NOISE formants
noteParamsGran = numNotesGran.collect({arg iGran;
	Array.with(startGran + (iGran * durGran), durGran+iGran, gain, risGran, decGran, noteFreqsGran.at(iGran))
});

// -------------------------------------------
// use -collect (iteration!) to pack into array, noteParams




// call noteParamsFunc to add notes to score!
noteParamsFunc.value(noteParams);
noteParamsFunc.value(noteParams1);
noteParamsFunc.value(noteParamsMid);
noteParamsFuncGran.value(synthDefGran, noteParamsGran, formDataGran);



// write score to sound file with the -write message
// NOTE: we're using argument names to specify the args. For 'duration', we're letting Ctk
//       do the work for us!
score.write(
    outputPath.standardizePath,
    sampleRate: s.sampleRate,
    headerFormat: headerFormat,
    sampleFormat: sampleFormat,
    options: ServerOptions.new.numOutputBusChannels_(numOutputChannels)
);
score.saveToFile("~/Desktop/Kedia_Final.scd".standardizePath);

)

SFPlayer("~/Desktop/Kedia_Final.wav".standardizePath).gui;
