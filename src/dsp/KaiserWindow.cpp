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

#include "KaiserWindow.h"

#include "MathUtilities.h"

KaiserWindow::Parameters
KaiserWindow::parametersForTransitionWidth(double attenuation,
					   double transition)
{
    Parameters p;
    p.length = 1 + (attenuation > 21.0 ?
		    ceil((attenuation - 7.95) / (2.285 * transition)) :
		    ceil(5.79 / transition));
    p.beta = (attenuation > 50.0 ? 
	      0.1102 * (attenuation - 8.7) :
	      attenuation > 21.0 ? 
	      0.5842 * pow(attenuation - 21.0, 0.4) + 0.07886 * (attenuation - 21.0) :
	      0);
    return p;
}

static double besselTerm(double x, int i)
{
    if (i == 0) {
	return 1;
    } else {
	double f = MathUtilities::factorial(i);
	return pow(x/2, i*2) / (f*f);
    }
}

static double bessel0(double x)
{
    double b = 0.0;
    for (int i = 0; i < 20; ++i) {
	b += besselTerm(x, i);
    }
    return b;
}

void
KaiserWindow::init()
{
    double denominator = bessel0(m_beta);
    bool even = (m_length % 2 == 0);
    for (int i = 0; i < (even ? m_length/2 : (m_length+1)/2); ++i) {
	double k = double(2*i) / double(m_length-1) - 1.0;
	m_window.push_back(bessel0(m_beta * sqrt(1.0 - k*k)) / denominator);
    }
    for (int i = 0; i < (even ? m_length/2 : (m_length-1)/2); ++i) {
        m_window.push_back(m_window[int(m_length/2) - i - 1]);
    }
}
