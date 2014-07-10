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

#ifndef CQ_KERNEL_H
#define CQ_KERNEL_H

#include "CQParameters.h"

#include <vector>
#include <complex>

class FFT;

class CQKernel
{
public:
    CQKernel(CQParameters params);
    ~CQKernel();

    bool isValid() const { return m_valid; }
    
    struct Properties {
        double sampleRate;
        double maxFrequency;
        double minFrequency;
        int binsPerOctave;
        int fftSize;
        int fftHop;
        int atomsPerFrame;
        int atomSpacing;
        int firstCentre;
        int lastCentre;
        double Q;
    };

    Properties getProperties() const { return m_p; }

    std::vector<std::complex<double> > processForward
        (const std::vector<std::complex<double> > &);

    std::vector<std::complex<double> > processInverse
        (const std::vector<std::complex<double> > &);

private:
    const CQParameters m_inparams;
    Properties m_p;
    bool m_valid;
    FFT *m_fft;

    struct KernelMatrix {
        std::vector<int> origin;
        std::vector<std::vector<std::complex<double> > > data;
    };
    KernelMatrix m_kernel;

    std::vector<double> makeWindow(int len) const;
    bool generateKernel();
    void finaliseKernel();
};

#endif
