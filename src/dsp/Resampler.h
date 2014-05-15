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

#ifndef RESAMPLER_H
#define RESAMPLER_H

#include <vector>

/**
 * Resampler resamples a stream from one integer sample rate to
 * another (arbitrary) rate, using a kaiser-windowed sinc filter.  The
 * results and performance are pretty similar to libraries such as
 * libsamplerate, though this implementation does not support
 * time-varying ratios (the ratio is fixed on construction).
 *
 * See also Decimator, which is faster and rougher but supports only
 * power-of-two downsampling factors.
 */
class Resampler
{
public:
    /**
     * Construct a Resampler to resample from sourceRate to
     * targetRate.
     */
    Resampler(int sourceRate, int targetRate);

    /**
     * Construct a Resampler to resample from sourceRate to
     * targetRate, using the given filter parameters.
     */
    Resampler(int sourceRate, int targetRate,
              double snr, double bandwidth);

    virtual ~Resampler();

    /**
     * Read n input samples from src and write resampled data to
     * dst. The return value is the number of samples written, which
     * will be no more than ceil((n * targetRate) / sourceRate). The
     * caller must ensure the dst buffer has enough space for the
     * samples returned.
     */
    int process(const double *src, double *dst, int n);

    /**
     * Read n input samples from src and return resampled data by
     * value.
     */
    std::vector<double> process(const double *src, int n);

    /**
     * Return the number of samples of latency at the output due by
     * the filter. (That is, the output will be delayed by this number
     * of samples relative to the input.)
     */
    int getLatency() const { return m_latency; }

    /**
     * Carry out a one-off resample of a single block of n
     * samples. The output is latency-compensated.
     */
    static std::vector<double> resample
    (int sourceRate, int targetRate, const double *data, int n);

private:
    int m_sourceRate;
    int m_targetRate;
    int m_gcd;
    int m_filterLength;
    int m_bufferLength;
    int m_latency;
    double m_peakToPole;
    
    struct Phase {
        int nextPhase;
        std::vector<double> filter;
        int drop;
    };

    Phase *m_phaseData;
    int m_phase;
    std::vector<double> m_buffer;
    int m_bufferOrigin;

    void initialise(double, double);
    double reconstructOne();
};

#endif
    
