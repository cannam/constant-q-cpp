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

#ifndef CQBASE_H
#define CQBASE_H

#include <vector>
#include <complex>

class CQBase // interface class
{
public:
    typedef std::complex<double> Complex;
    typedef std::vector<double> RealSequence;
    typedef std::vector<double> RealColumn;
    typedef std::vector<Complex> ComplexSequence;
    typedef std::vector<Complex> ComplexColumn;
    typedef std::vector<RealColumn> RealBlock;
    typedef std::vector<ComplexColumn> ComplexBlock;

    virtual double getSampleRate() const = 0;
    virtual int getBinsPerOctave() const = 0;
    virtual int getOctaves() const = 0; 
    virtual int getTotalBins() const = 0;
    virtual int getColumnHop() const = 0;
    virtual int getLatency() const = 0;
    virtual double getMaxFrequency() const = 0;
    virtual double getMinFrequency() const = 0; // actual min, not that passed to ctor
    virtual double getBinFrequency(int bin) const = 0;
};

#endif
