/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/*
    Constant-Q library
    Copyright (c) 2013-2014 Queen Mary, University of London

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    Except as contained in this notice, the names of the Centre for
    Digital Music; Queen Mary, University of London; and Chris Cannam
    shall not be used in advertising or otherwise to promote the sale,
    use or other dealings in this Software without prior written
    authorization.
*/

#ifndef CQ_PARAMETERS_H
#define CQ_PARAMETERS_H

/**
 * Common parameters for constructing Constant-Q implementation
 * objects (both forward and inverse transforms).
 */
class CQParameters
{
public:
    enum WindowType {
	SqrtBlackmanHarris,
	SqrtBlackman,
	SqrtHann,
	BlackmanHarris,
	Blackman,
	Hann,
    };

    enum DecimatorType {
        BetterDecimator,
        FasterDecimator
    };
    
    /**
     * Construct a set of parameters with the given input signal
     * sample rate, frequency range, and number of bins per
     * octave. The remaining parameters will take their usual
     * defaults; if you want to change them, just assign the
     * respective data members after construction.
     */
    CQParameters(double _sampleRate, 
		 double _minFrequency, 
		 double _maxFrequency,
		 int _binsPerOctave) :
	sampleRate(_sampleRate),
	minFrequency(_minFrequency),
	maxFrequency(_maxFrequency),
	binsPerOctave(_binsPerOctave),
	q(1.0),                     // Q scaling factor
	atomHopFactor(0.25),        // hop size of shortest temporal atom
	threshold(0.0005),          // sparsity threshold for resulting kernel
	window(SqrtBlackmanHarris), // window shape
        decimator(BetterDecimator)  // decimator quality setting
    { }

    /**
     * Sampling rate of input signal.
     */
    double sampleRate;

    /**
     * Minimum frequency desired to include in Constant-Q output. The
     * actual minimum will normally be calculated as a round number of
     * octaves below the maximum frequency, and may differ from this.
     */
    double minFrequency;

    /**
     * Maximum frequency to include in Constant-Q output.
     */
    double maxFrequency;

    /**
     * Number of output frequency bins per octave.
     */
    int binsPerOctave;

    /**
     * Spectral atom bandwidth scaling factor. q == 1 is optimal for
     * reconstruction, q < 1 increases redundancy (smearing) in the
     * frequency domain but improves time resolution.
     */
    double q;

    /**
     * Hop size between temporal atoms, where 1 == no overlap and
     * smaller values indicate overlapping atoms.
     */
    double atomHopFactor;

    /**
     * Sparsity threshold for Constant-Q kernel: values with magnitude
     * smaller than this are truncated to zero.
     */
    double threshold;

    /**
     * Window shape to use for the Constant-Q kernel atoms.
     */
    WindowType window;

    /**
     * Quality setting for the sample rate decimator.
     */
    DecimatorType decimator;
};

#endif

