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

#ifndef WINDOW_H
#define WINDOW_H

#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#include "pi.h"

enum WindowType {
    RectangularWindow,
    BartlettWindow,
    HammingWindow,
    HanningWindow,
    BlackmanWindow,
    BlackmanHarrisWindow,

    FirstWindow = RectangularWindow,
    LastWindow = BlackmanHarrisWindow
};

/**
 * Various shaped windows for sample frame conditioning, including
 * cosine windows (Hann etc) and triangular and rectangular windows.
 */
template <typename T>
class Window
{
public:
    /**
     * Construct a windower of the given type and size. 
     *
     * Note that the cosine windows are periodic by design, rather
     * than symmetrical. (A window of size N is equivalent to a
     * symmetrical window of size N+1 with the final element missing.)
     */
    Window(WindowType type, int size) : m_type(type), m_size(size) { encache(); }
    Window(const Window &w) : m_type(w.m_type), m_size(w.m_size) { encache(); }
    Window &operator=(const Window &w) {
	if (&w == this) return *this;
	m_type = w.m_type;
	m_size = w.m_size;
	encache();
	return *this;
    }
    virtual ~Window() { delete[] m_cache; }
    
    void cut(T *src) const { cut(src, src); }
    void cut(const T *src, T *dst) const {
	for (int i = 0; i < m_size; ++i) dst[i] = src[i] * m_cache[i];
    }

    WindowType getType() const { return m_type; }
    int getSize() const { return m_size; }

    std::vector<T> getWindowData() const {
        std::vector<T> d;
        for (int i = 0; i < m_size; ++i) {
            d.push_back(m_cache[i]);
        }
        return d;
    }

protected:
    WindowType m_type;
    int m_size;
    T *m_cache;
    
    void encache();
};

template <typename T>
void Window<T>::encache()
{
    int n = m_size;
    T *mult = new T[n];
    int i;
    for (i = 0; i < n; ++i) mult[i] = 1.0;

    switch (m_type) {
		
    case RectangularWindow:
        for (i = 0; i < n; ++i) {
            mult[i] = mult[i] * 0.5;
	}
	break;
	    
    case BartlettWindow:
        if (n == 2) {
            mult[0] = mult[1] = 0; // "matlab compatible"
        } else if (n == 3) {
            mult[0] = 0;
            mult[1] = mult[2] = 2./3.;
        } else if (n > 3) {
            for (i = 0; i < n/2; ++i) {
                mult[i] = mult[i] * (i / T(n/2));
                mult[i + n - n/2] = mult[i + n - n/2] * (1.0 - (i / T(n/2)));
            }
	}
	break;
	    
    case HammingWindow:
        if (n > 1) {
            for (i = 0; i < n; ++i) {
                mult[i] = mult[i] * (0.54 - 0.46 * cos(2 * M_PI * i / n));
            }
	}
	break;
	    
    case HanningWindow:
        if (n > 1) {
            for (i = 0; i < n; ++i) {
                mult[i] = mult[i] * (0.50 - 0.50 * cos(2 * M_PI * i / n));
            }
	}
	break;
	    
    case BlackmanWindow:
        if (n > 1) {
            for (i = 0; i < n; ++i) {
                mult[i] = mult[i] * (0.42 - 0.50 * cos(2 * M_PI * i / n)
                                     + 0.08 * cos(4 * M_PI * i / n));
            }
	}
	break;
	    
    case BlackmanHarrisWindow:
        if (n > 1) {
            for (i = 0; i < n; ++i) {
                mult[i] = mult[i] * (0.35875
                                     - 0.48829 * cos(2 * M_PI * i / n)
                                     + 0.14128 * cos(4 * M_PI * i / n)
                                     - 0.01168 * cos(6 * M_PI * i / n));
            }
	}
	break;
    }
	   
    m_cache = mult;
}

#endif
