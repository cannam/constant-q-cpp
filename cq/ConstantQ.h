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

#ifndef CONSTANTQ_H
#define CONSTANTQ_H

#include "CQBase.h"
#include "CQParameters.h"
#include "CQKernel.h"

class Resampler;
class FFTReal;

/**
 * Calculate a complex sparse constant-Q representation from
 * time-domain input. The input of each \ref process call is a single
 * frame of time-domain samples; the output is a series of columns of
 * varying height. See \ref process for details.
 *
 * For a real (magnitude-only) interpolated dense representation, see
 * CQSpectrogram.
 */
class ConstantQ : public CQBase
{
public:
    /**
     * Construct a complex Constant-Q transform object using the given
     * transform parameters.
     */
    ConstantQ(CQParameters params);
    virtual ~ConstantQ();

    // CQBase methods, see CQBase.h for documentation
    virtual bool isValid() const { return m_kernel && m_kernel->isValid(); }
    virtual double getSampleRate() const { return m_sampleRate; }
    virtual int getBinsPerOctave() const { return m_binsPerOctave; }
    virtual int getOctaves() const { return m_octaves; }
    virtual int getTotalBins() const { return m_octaves * m_binsPerOctave; }
    virtual int getColumnHop() const { return m_p.fftHop / m_p.atomsPerFrame; }
    virtual int getLatency() const { return m_outputLatency; } 
    virtual double getMaxFrequency() const { return m_p.maxFrequency; }
    virtual double getMinFrequency() const;
    virtual double getBinFrequency(double bin) const; // bin may be nonintegral

    /**
     * Given a series of time-domain samples, return a series of
     * constant-Q columns. Any samples left over (that did not fit
     * into a constant-Q processing block) are saved for the next call
     * to process or getRemainingBlocks. 
     *
     * The input is assumed to be a single frame of time-domain sample
     * values, such that consecutive calls to \ref process receive
     * contiguous frames from the source signal. Each frame may be of
     * any length in samples.
     *
     * Each output column contains a series of constant-Q bin values
     * ordered from highest to lowest frequency.
     *
     * Output columns are of varying height: each will contain at
     * least getBinsPerOctave() values, because the highest-frequency
     * octave is always present, but a second octave (if requested)
     * will appear only in alternate columns, a third octave only in
     * every fourth column, and so on.
     *
     * If you need a format in which all columns are of equal height
     * and every bin contains a value, use \ref CQSpectrogram instead
     * of ConstantQ.
     */
    ComplexBlock process(const RealSequence &);

    /**
     * Return the remaining constant-Q columns following the end of
     * processing. Any buffered input is padded so as to ensure that
     * all input provided to process() will have been returned.
     */
    ComplexBlock getRemainingOutput();

private:
    const CQParameters m_inparams;
    const double m_sampleRate;
    const double m_maxFrequency;
    const double m_minFrequency;
    const int m_binsPerOctave;

    int m_octaves;
    CQKernel *m_kernel;
    CQKernel::Properties m_p;
    int m_bigBlockSize;

    std::vector<Resampler *> m_decimators;
    std::vector<RealSequence> m_buffers;

    int m_outputLatency;

    FFTReal *m_fft;

    void initialise();
    ComplexBlock processOctaveBlock(int octave);
};

#endif

