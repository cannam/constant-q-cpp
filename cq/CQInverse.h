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

#ifndef CQINVERSE_H
#define CQINVERSE_H

#include "CQBase.h"
#include "CQKernel.h"

class Resampler;
class FFTReal;

/**
 * Calculate an inverse constant-Q transform. The input must be the
 * same representation as returned as output of a \ref ConstantQ
 * object with the same parameters. The output is a time-domain
 * signal.
 *
 * Note that you cannot perform an inverse transform from the
 * magnitude-only output of \ref CQSpectrogram; you need the complex
 * valued data from \ref ConstantQ.
 *
 * Our implementation of the Constant-Q transform is not exactly
 * invertible, and this produces only an approximation of the original
 * signal (see publications for details).
 */
class CQInverse : public CQBase
{
public:
    /**
     * Construct an inverse Constant-Q transform object using the
     * given transform parameters.
     */
    CQInverse(CQParameters params);
    virtual ~CQInverse();

    // CQBase methods, see CQBase.h for documentation
    virtual bool isValid() const { return m_kernel && m_kernel->isValid(); }
    virtual double getSampleRate() const { return m_sampleRate; }
    virtual int getBinsPerOctave() const { return m_binsPerOctave; }
    virtual int getOctaves() const { return m_octaves; }
    virtual int getTotalBins() const { return m_octaves * m_binsPerOctave; }
    virtual int getColumnHop() const { return m_p.fftHop / m_p.atomsPerFrame; }
    virtual int getLatency() const { return m_outputLatency; } 
    virtual double getMaxFrequency() const { return m_p.maxFrequency; }
    virtual double getMinFrequency() const; // actual min, not that passed to ctor
    virtual double getBinFrequency(double bin) const;

    /**
     * Given a series of constant-Q columns in the form produced by
     * the \ref ConstantQ class, return a series of time-domain
     * samples resulting from approximately inverting the constant-Q
     * transform.
     */
    RealSequence process(const ComplexBlock &);

    /**
     * Return the remaining time-domain samples following the end of
     * processing.
     */
    RealSequence getRemainingOutput();

private:
    const CQParameters m_inparams;
    const double m_sampleRate;
    const double m_maxFrequency;
    const double m_minFrequency;
    const int m_binsPerOctave;

    int m_octaves;
    CQKernel *m_kernel;
    CQKernel::Properties m_p;

    std::vector<Resampler *> m_upsamplers;
    std::vector<RealSequence> m_buffers;
    std::vector<RealSequence> m_olaBufs; // fixed-length, for overlap-add
    
    int m_outputLatency;

    FFTReal *m_fft;
    
    void initialise();
    void processOctave(int octave, const ComplexBlock &block);
    void processOctaveColumn(int octave, const ComplexColumn &column);
    void overlapAddAndResample(int octave, const RealSequence &);
    RealSequence drawFromBuffers();
};

#endif
