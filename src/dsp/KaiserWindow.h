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

#ifndef KAISER_WINDOW_H
#define KAISER_WINDOW_H

#include <vector>
#include <cmath>

#include "pi.h"

/**
 * Kaiser window: A windower whose bandwidth and sidelobe height
 * (signal-noise ratio) can be specified. These parameters are traded
 * off against the window length.
 */
class KaiserWindow
{
public:
    struct Parameters {
	int length;
	double beta;
    };

    /**
     * Construct a Kaiser windower with the given length and beta
     * parameter.
     */
    KaiserWindow(Parameters p) : m_length(p.length), m_beta(p.beta) { init(); }

    /**
     * Construct a Kaiser windower with the given attenuation in dB
     * and transition width in samples.
     */
    static KaiserWindow byTransitionWidth(double attenuation,
					  double transition) {
	return KaiserWindow
	    (parametersForTransitionWidth(attenuation, transition));
    }

    /**
     * Construct a Kaiser windower with the given attenuation in dB
     * and transition bandwidth in Hz for the given samplerate.
     */
    static KaiserWindow byBandwidth(double attenuation,
				    double bandwidth,
				    double samplerate) {
	return KaiserWindow
	    (parametersForBandwidth(attenuation, bandwidth, samplerate));
    }

    /**
     * Obtain the parameters necessary for a Kaiser window of the
     * given attenuation in dB and transition width in samples.
     */
    static Parameters parametersForTransitionWidth(double attenuation,
						   double transition);

    /**
     * Obtain the parameters necessary for a Kaiser window of the
     * given attenuation in dB and transition bandwidth in Hz for the
     * given samplerate.
     */
    static Parameters parametersForBandwidth(double attenuation,
					     double bandwidth,
					     double samplerate) {
	return parametersForTransitionWidth
	    (attenuation, (bandwidth * 2 * M_PI) / samplerate);
    } 

    int getLength() const {
	return m_length;
    }

    const double *getWindow() const { 
	return m_window.data();
    }

    void cut(double *src) const { 
	cut(src, src); 
    }

    void cut(const double *src, double *dst) const {
	for (int i = 0; i < m_length; ++i) {
	    dst[i] = src[i] * m_window[i];
	}
    }

private:
    int m_length;
    double m_beta;
    std::vector<double> m_window;

    void init();
};

#endif
