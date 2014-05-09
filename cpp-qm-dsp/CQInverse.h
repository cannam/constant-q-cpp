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

#ifndef CQ_INVERSE_H
#define CQ_INVERSE_H

#include "CQKernel.h"

#include <vector>

class Resampler;
class FFTReal;

class CQInverse
{
public:
    CQInverse(double sampleRate,
	      double minFreq, double maxFreq,
	      int binsPerOctave);
    ~CQInverse();

    double getSampleRate() const { return m_sampleRate; }
    int getBinsPerOctave() const { return m_binsPerOctave; }
    int getOctaves() const { return m_octaves; }
    int getTotalBins() const { return m_octaves * m_binsPerOctave; }
    int getColumnHop() const { return m_p.fftHop / m_p.atomsPerFrame; }
    int getLatency() const { return m_outputLatency; } 
    double getMaxFrequency() const { return m_p.maxFrequency; }
    double getMinFrequency() const; // actual min, not that passed to ctor
    double getBinFrequency(int bin) const;

    // Input is the format produced by ConstantQ class, not
    // CQInterpolated (or can we make this support either?)

    //!!! no, we need complex not magnitudes! CQ should probably
    //!!! produce totally raw output and something like CQInterpolated
    //!!! do the magnitude stuff as well as interpolation
    std::vector<double> process(const std::vector<std::vector<double> > &);
    std::vector<double> getRemainingOutput();

private:
    double m_sampleRate;
    double m_maxFrequency;
    double m_minFrequency;
    int m_binsPerOctave;
    int m_octaves;

    CQKernel *m_kernel;
    CQKernel::Properties m_p;
    int m_bigBlockSize;

    std::vector<Resampler *> m_upsamplers;
    std::vector<std::vector<double> > m_buffers;
    
    int m_outputLatency;

    FFTReal *m_fft;
    
    void initialise();
};

#endif
