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

#include "CQInverse.h"

#include "CQKernel.h"

#include "dsp/rateconversion/Resampler.h"
#include "maths/MathUtilities.h"
#include "dsp/transforms/FFT.h"

using std::vector;
using std::complex;
using std::cerr;
using std::endl;

typedef std::complex<double> C;

CQInverse::CQInverse(double sampleRate,
                     double minFreq,
                     double maxFreq,
                     int binsPerOctave) :
    m_sampleRate(sampleRate),
    m_maxFrequency(maxFreq),
    m_minFrequency(minFreq),
    m_binsPerOctave(binsPerOctave),
    m_fft(0)
{
    if (minFreq <= 0.0 || maxFreq <= 0.0) {
        throw std::invalid_argument("Frequency extents must be positive");
    }

    initialise();
}

CQInverse::~CQInverse()
{
    delete m_fft;
    for (int i = 0; i < (int)m_upsamplers.size(); ++i) {
        delete m_upsamplers[i];
    }
    delete m_kernel;
}

double
CQInverse::getMinFrequency() const
{
    return m_p.minFrequency / pow(2.0, m_octaves - 1);
}

double
CQInverse::getBinFrequency(int bin) const
{
    return getMinFrequency() * pow(2, (double(bin) / getBinsPerOctave()));
}

void
CQInverse::initialise()
{
    m_octaves = int(ceil(log2(m_maxFrequency / m_minFrequency)));
    m_kernel = new CQKernel(m_sampleRate, m_maxFrequency, m_binsPerOctave);
    m_p = m_kernel->getProperties();
    
    // Use exact powers of two for resampling rates. They don't have
    // to be related to our actual samplerate: the resampler only
    // cares about the ratio, but it only accepts integer source and
    // target rates, and if we start from the actual samplerate we
    // risk getting non-integer rates for lower octaves

    int sourceRate = pow(2, m_octaves);
    vector<int> latencies;

    // top octave, no resampling
    latencies.push_back(0);
    m_upsamplers.push_back(0);

    for (int i = 1; i < m_octaves; ++i) {

        int factor = pow(2, i);

        Resampler *r = new Resampler
            (sourceRate / factor, sourceRate, 60, 0.02);

	// See ConstantQ.cpp for discussion on latency -- output
	// latency here is at target rate which, this way around, is
	// what we want

        latencies.push_back(r->getLatency());
        m_upsamplers.push_back(r);
    }

    m_bigBlockSize = m_p.fftSize * pow(2, m_octaves - 1);


//!!! GOT HERE!
#error continue from here please!
/*
    // Now add in the extra padding and compensate for hops that must
    // be dropped in order to align the atom centres across
    // octaves. Again this is a bit trickier because we are doing it
    // at input rather than output and so must work in per-octave
    // sample rates rather than output blocks

    int emptyHops = m_p.firstCentre / m_p.atomSpacing;

    vector<int> drops;
    for (int i = 0; i < m_octaves; ++i) {
	int factor = pow(2, i);
	int dropHops = emptyHops * pow(2, m_octaves - i - 1) - emptyHops;
	int drop = ((dropHops * m_p.fftHop) * factor) / m_p.atomsPerFrame;
	drops.push_back(drop);
    }

    int maxLatPlusDrop = 0;
    for (int i = 0; i < m_octaves; ++i) {
	int latPlusDrop = latencies[i] + drops[i];
	if (latPlusDrop > maxLatPlusDrop) maxLatPlusDrop = latPlusDrop;
    }

    // we want to design totalLatency such that totalLatency -
    // latencies[0] - drops[0] is a multiple of m_p.fftHop, so that we
    // can get identical results in octave 0 to our reference
    // implementation, making for easier testing (though other octaves
    // will differ because of different resampler implementations)

    int totalLatency = maxLatPlusDrop;
    int lat0 = totalLatency - latencies[0] - drops[0];
    totalLatency = ceil(double(lat0 / m_p.fftHop) * m_p.fftHop)
	+ latencies[0] + drops[0];

//    cerr << "total latency = " << totalLatency << endl;

    // Padding as in the reference (will be introduced with the
    // latency compensation in the loop below)
    m_outputLatency = totalLatency + m_bigBlockSize
	- m_p.firstCentre * pow(2, m_octaves-1);

//    cerr << "m_bigBlockSize = " << m_bigBlockSize << ", firstCentre = "
//	 << m_p.firstCentre << ", m_octaves = " << m_octaves << ", so m_outputLatency = " << m_outputLatency << endl;

    for (int i = 0; i < m_octaves; ++i) {

	double factor = pow(2, i);

	// Calculate the difference between the total latency applied
	// across all octaves, and the existing latency due to the
	// decimator for this octave, and then convert it back into
	// the sample rate appropriate for the output latency of this
	// decimator -- including one additional big block of padding
	// (as in the reference).

	double octaveLatency =
	    double(totalLatency - latencies[i] - drops[i]
		   + m_bigBlockSize) / factor;

        m_buffers.push_back
            (vector<double>(int(round(octaveLatency)), 0.0));
    }
*/
    m_fft = new FFTReal(m_p.fftSize);
}
