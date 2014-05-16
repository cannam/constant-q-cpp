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

    CQParameters(double _sampleRate, 
		 double _minFrequency, 
		 double _maxFrequency,
		 int _binsPerOctave) :
	sampleRate(_sampleRate),
	minFrequency(_minFrequency),
	maxFrequency(_maxFrequency),
	binsPerOctave(_binsPerOctave),
	q(1.0),                    // Q scaling factor
	atomHopFactor(0.25),       // hop size of shortest temporal atom
	threshold(0.0005),         // sparsity threshold for resulting kernel
	window(SqrtBlackmanHarris) // window shape
    { }

    double sampleRate;
    double minFrequency;
    double maxFrequency;
    int binsPerOctave;

    double q;
    double atomHopFactor;
    double threshold;
    WindowType window;
};

#endif


    
