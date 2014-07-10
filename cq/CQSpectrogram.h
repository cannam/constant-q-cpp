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

#ifndef CQSPECTROGRAM_H
#define CQSPECTROGRAM_H

#include "ConstantQ.h"

/**
 * Calculate a dense constant-Q magnitude spectrogram from time-domain
 * input. The input of each \ref process call is a single frame of
 * time-domain samples; the output is a series of fixed-height
 * columns. See \ref process for details.
 *
 * If you need the full complex-valued constant-Q output, you must use
 * the \ref ConstantQ class instead.
 */
class CQSpectrogram : public CQBase
{
public:
    enum Interpolation {
        /// leave empty cells as zero
	InterpolateZeros,
        /// replace empty cells with a repeat of the previous column
	InterpolateHold,
        /// perform linear interpolation between consecutive time cells
	InterpolateLinear,
    };

    /**
     * Construct a Constant-Q magnitude spectrogram object using the
     * given transform parameters.
     */
    CQSpectrogram(CQParameters params, Interpolation interpolation);
    virtual ~CQSpectrogram();

    // CQBase methods, see CQBase.h for documentation
    virtual bool isValid() const { return m_cq.isValid(); }
    virtual double getSampleRate() const { return m_cq.getSampleRate(); }
    virtual int getBinsPerOctave() const { return m_cq.getBinsPerOctave(); }
    virtual int getOctaves() const { return m_cq.getOctaves(); }
    virtual int getTotalBins() const { return m_cq.getTotalBins(); }
    virtual int getColumnHop() const { return m_cq.getColumnHop(); }
    virtual int getLatency() const { return m_cq.getLatency(); } 
    virtual double getMaxFrequency() const { return m_cq.getMaxFrequency(); }
    virtual double getMinFrequency() const { return m_cq.getMinFrequency(); }
    virtual double getBinFrequency(double bin) const { return m_cq.getBinFrequency(bin); }

    /**
     * Given a series of time-domain samples, return a series of
     * constant-Q magnitude columns. Any samples left over (that did
     * not fit into a constant-Q processing block) are saved for the
     * next call to process or getRemainingBlocks.
     *
     * The input is assumed to be a single frame of time-domain sample
     * values, such that consecutive calls to \ref process receive
     * contiguous frames from the source signal. Each frame may be of
     * any length in samples.
     *
     * Each output column contains a series of constant-Q bin value
     * magnitudes, ordered from highest to lowest frequency.
     *  
     * The columns are all of the same height, but they might not all
     * be populated, depending on the interpolation mode: in
     * InterpolateZeros mode, the lower octaves (which are spaced more
     * widely in the raw constant-Q than the highest octave) will
     * contain zeros for the undefined values, but in the other
     * interpolation modes every cell will be filled.
     *
     * To obtain raw, complex constant-Q bin values, use the ConstantQ
     * class.
     */
    RealBlock process(const RealSequence &);

    /**
     * Return the remaining constant-Q magnitude columns following the
     * end of processing. Any buffered input is padded so as to ensure
     * that all input provided to process() will have been returned.
     */
    RealBlock getRemainingOutput();

private:
    ConstantQ m_cq;
    Interpolation m_interpolation;

    RealBlock m_buffer;
    RealBlock postProcess(const ComplexBlock &, bool insist);
    RealBlock fetchHold(bool insist);
    RealBlock fetchLinear(bool insist);
    RealBlock linearInterpolated(const RealBlock &, int, int);
    RealColumn m_prevColumn;
};

#endif
