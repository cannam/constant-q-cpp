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

#include "SincWindow.h"

#include <cmath>

void
SincWindow::init()
{
    if (m_length < 1) {
	return;
    } else if (m_length < 2) {
	m_window.push_back(1);
	return;
    } else {

	int n0 = (m_length % 2 == 0 ? m_length/2 : (m_length - 1)/2);
	int n1 = (m_length % 2 == 0 ? m_length/2 : (m_length + 1)/2);
	double m = 2 * M_PI / m_p;

	for (int i = 0; i < n0; ++i) {
	    double x = ((m_length / 2) - i) * m;
	    m_window.push_back(sin(x) / x);
	}

	m_window.push_back(1.0);

	for (int i = 1; i < n1; ++i) {
	    double x = i * m;
	    m_window.push_back(sin(x) / x);
	}
    }
}

