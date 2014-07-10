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

/**
 * Interface class for Constant-Q implementations, containing common
 * type declarations and means to query configuration parameters.
 */
class CQBase
{
public:
    /// A single complex-valued sample.
    typedef std::complex<double> Complex;

    /// A series of real-valued samples ordered in time.
    typedef std::vector<double> RealSequence;

    /// A series of real-valued samples ordered by bin (frequency or similar).
    typedef std::vector<double> RealColumn;

    /// A series of complex-valued samples ordered in time.
    typedef std::vector<Complex> ComplexSequence;

    /// A series of complex-valued samples ordered by bin (frequency or similar).
    typedef std::vector<Complex> ComplexColumn;

    /// A matrix of real-valued samples, indexed by time then bin number.
    typedef std::vector<RealColumn> RealBlock;

    /// A matrix of complex-valued samples, indexed by time then bin number.
    typedef std::vector<ComplexColumn> ComplexBlock;

    /**
     * Return true if the Constant-Q implementation was successfully
     * constructed, with a valid set of initialisation parameters.
     */
    virtual bool isValid() const = 0;
    
    /**
     * Return the sample rate used when constructing the specific
     * Constant-Q implementation.
     */
    virtual double getSampleRate() const = 0;

    /**
     * Return the number of bins per octave specified when
     * constructing the Constant-Q implementation.
     */
    virtual int getBinsPerOctave() const = 0;

    /**
     * Return the number of octaves spanned by the Constant-Q
     * transform.
     */
    virtual int getOctaves() const = 0; 

    /**
     * Return the total number of bins in each Constant-Q column
     * (i.e. bins-per-octave times number of octaves).
     */
    virtual int getTotalBins() const = 0;

    /**
     * Return the spacing, in samples at the sample rate returned from
     * getSampleRate(), between one column and the next.
     */
    virtual int getColumnHop() const = 0;

    /**
     * Return the latency of Constant-Q calculation, in samples at the
     * sample rate returned from getSampleRate().
     */
    virtual int getLatency() const = 0;

    /**
     * Return the maximum frequency of the Constant-Q output, i.e. the
     * frequency of the highest bin in the output. This will normally
     * be the same as the maximum frequency passed to the constructor
     * of the specific Constant-Q implementation.
     */
    virtual double getMaxFrequency() const = 0;

    /**
     * Return the minimum frequency of the Constant-Q output, i.e. the
     * frequency of the lowest bin in the output. This is derived from
     * the maximum frequency and octave count, and is not necessarily
     * the same as any minimum frequency requested when constructing
     * the Constant-Q implementation.
     */
    virtual double getMinFrequency() const = 0;

    /**
     * Return the frequency of a given bin in the Constant-Q
     * output. This actually maps a continuous "bin scale" value to 
     * frequency: the bin parameter does not have to be an integer.
     */
    virtual double getBinFrequency(double bin) const = 0;
};

#endif
