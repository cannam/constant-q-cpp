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

#ifndef CQ_INTERPOLATED_H
#define CQ_INTERPOLATED_H

#include <vector>
#include "ConstantQ.h"

class CQInterpolated
{
public:
    enum Interpolation {
	None, // leave empty cells empty
	Hold, // repeat prior cell
	Linear, // linear interpolation between consecutive time cells
    };

    CQInterpolated(double sampleRate,
		   double minFreq, double maxFreq,
		   int binsPerOctave,
		   Interpolation interpolation);
    ~CQInterpolated();

    double getSampleRate() const { return m_cq.getSampleRate(); }
    int getBinsPerOctave() const { return m_cq.getBinsPerOctave(); }
    int getOctaves() const { return m_cq.getOctaves(); }
    int getTotalBins() const { return m_cq.getTotalBins(); }
    int getColumnHop() const { return m_cq.getColumnHop(); }
    int getLatency() const { return m_cq.getLatency(); } 
    double getMaxFrequency() const { return m_cq.getMaxFrequency(); }
    double getMinFrequency() const { return m_cq.getMinFrequency(); }
    double getBinFrequency(int bin) const { return m_cq.getBinFrequency(bin); }

    std::vector<std::vector<double> > process(const std::vector<double> &);
    std::vector<std::vector<double> > getRemainingBlocks();

private:
    ConstantQ m_cq;
    Interpolation m_interpolation;

    typedef std::vector<std::vector<double> > Grid;
    Grid m_buffer;
    Grid postProcess(Grid, bool insist);
    Grid fetchHold(bool insist);
    Grid fetchLinear(bool insist);
    std::vector<double> m_prevColumn;
};

#endif
