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

class CQSpectrogram : public CQBase
{
public:
    enum Interpolation {
	InterpolateZeros, // leave empty cells as zero
	InterpolateHold, // repeat prior cell
	InterpolateLinear, // linear interpolation between consecutive time cells
    };

    CQSpectrogram(double sampleRate,
                  double minFreq, double maxFreq,
                  int binsPerOctave,
                  Interpolation interpolation);
    ~CQSpectrogram();

    virtual double getSampleRate() const { return m_cq.getSampleRate(); }
    virtual int getBinsPerOctave() const { return m_cq.getBinsPerOctave(); }
    virtual int getOctaves() const { return m_cq.getOctaves(); }
    virtual int getTotalBins() const { return m_cq.getTotalBins(); }
    virtual int getColumnHop() const { return m_cq.getColumnHop(); }
    virtual int getLatency() const { return m_cq.getLatency(); } 
    virtual double getMaxFrequency() const { return m_cq.getMaxFrequency(); }
    virtual double getMinFrequency() const { return m_cq.getMinFrequency(); }
    virtual double getBinFrequency(int bin) const { return m_cq.getBinFrequency(bin); }

    RealBlock process(const RealSequence &);
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
