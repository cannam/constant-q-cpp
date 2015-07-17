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

#include "ConstantQ.h"

#include "CQKernel.h"

#include "dsp/Resampler.h"
#include "dsp/MathUtilities.h"
#include "dsp/FFT.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <cmath>

using std::vector;
using std::cerr;
using std::endl;

//#define DEBUG_CQ 1

ConstantQ::ConstantQ(CQParameters params) :
    m_inparams(params),
    m_sampleRate(params.sampleRate),
    m_maxFrequency(params.maxFrequency),
    m_minFrequency(params.minFrequency),
    m_binsPerOctave(params.binsPerOctave),
    m_kernel(0),
    m_fft(0)
{
    if (m_minFrequency <= 0.0 || m_maxFrequency <= 0.0) {
        throw std::invalid_argument("Frequency extents must be positive");
    }

    initialise();
}

ConstantQ::~ConstantQ()
{
    delete m_fft;
    for (int i = 0; i < (int)m_decimators.size(); ++i) {
        delete m_decimators[i];
    }
    delete m_kernel;
}

double
ConstantQ::getMinFrequency() const
{
    return m_p.minFrequency / pow(2.0, m_octaves - 1);
}

double
ConstantQ::getBinFrequency(double bin) const
{
    // our bins are returned in high->low order
    bin = (getBinsPerOctave() * getOctaves()) - bin - 1;
    return getMinFrequency() * pow(2, (bin / getBinsPerOctave()));
}

void
ConstantQ::initialise()
{
    m_octaves = int(ceil(log(m_maxFrequency / m_minFrequency) / log(2)));

    if (m_octaves < 1) {
        m_kernel = 0; // incidentally causing isValid() to return false
        return;
    }

    m_kernel = new CQKernel(m_inparams);
    m_p = m_kernel->getProperties();
    
    if (!m_kernel->isValid()) {
        return;
    }

    // Use exact powers of two for resampling rates. They don't have
    // to be related to our actual samplerate: the resampler only
    // cares about the ratio, but it only accepts integer source and
    // target rates, and if we start from the actual samplerate we
    // risk getting non-integer rates for lower octaves

    int sourceRate = pow(2, m_octaves);
    vector<int> latencies;

    // top octave, no resampling
    latencies.push_back(0);
    m_decimators.push_back(0);

    for (int i = 1; i < m_octaves; ++i) {

        int factor = pow(2, i);

        Resampler *r;

        if (m_inparams.decimator == CQParameters::BetterDecimator) {
            r = new Resampler
                (sourceRate, sourceRate / factor, 50, 0.05);
        } else {
            r = new Resampler
                (sourceRate, sourceRate / factor, 25, 0.3);
        }                

#ifdef DEBUG_CQ
        cerr << "forward: octave " << i << ": resample from " << sourceRate << " to " << sourceRate / factor << endl;
#endif

        // We need to adapt the latencies so as to get the first input
        // sample to be aligned, in time, at the decimator output
        // across all octaves.
        // 
        // Our decimator uses a linear phase filter, but being causal
        // it is not zero phase: it has a latency that depends on the
        // decimation factor. Those latencies have been calculated
        // per-octave and are available to us in the latencies
        // array. Left to its own devices, the first input sample will
        // appear at output sample 0 in the highest octave (where no
        // decimation is needed), sample number latencies[1] in the
        // next octave down, latencies[2] in the next one, etc. We get
        // to apply some artificial per-octave latency after the
        // decimator in the processing chain, in order to compensate
        // for the differing latencies associated with different
        // decimation factors. How much should we insert?
        //
        // The outputs of the decimators are at different rates (in
        // terms of the relation between clock time and samples) and
        // we want them aligned in terms of time. So, for example, a
        // latency of 10 samples with a decimation factor of 2 is
        // equivalent to a latency of 20 with no decimation -- they
        // both result in the first output sample happening at the
        // same equivalent time in milliseconds.
	// 
	// So here we record the latency added by the decimator, in
	// terms of the sample rate of the undecimated signal. Then we
	// use that to compensate in a moment, when we've discovered
	// what the longest latency across all octaves is.

        latencies.push_back(r->getLatency() * factor);
        m_decimators.push_back(r);
    }

    m_bigBlockSize = m_p.fftSize * pow(2, m_octaves - 1);

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

    int totalLatency = maxLatPlusDrop;

    int lat0 = totalLatency - latencies[0] - drops[0];
    totalLatency = ceil(double(lat0 / m_p.fftHop) * m_p.fftHop)
	+ latencies[0] + drops[0];

    // We want (totalLatency - latencies[i]) to be a multiple of 2^i
    // for each octave i, so that we do not end up with fractional
    // octave latencies below. In theory this is hard, in practice if
    // we ensure it for the last octave we should be OK.
    double finalOctLat = latencies[m_octaves-1];
    double finalOctFact = pow(2, m_octaves-1);
    totalLatency =
        int(finalOctLat +
            finalOctFact *
            ceil((totalLatency - finalOctLat) / finalOctFact) + .5);

#ifdef DEBUG_CQ
    cerr << "total latency = " << totalLatency << endl;
#endif

    // Padding as in the reference (will be introduced with the
    // latency compensation in the loop below)
    m_outputLatency = totalLatency + m_bigBlockSize
	- m_p.firstCentre * pow(2, m_octaves-1);

#ifdef DEBUG_CQ
    cerr << "m_bigBlockSize = " << m_bigBlockSize << ", firstCentre = "
	 << m_p.firstCentre << ", m_octaves = " << m_octaves
         << ", so m_outputLatency = " << m_outputLatency << endl;
#endif

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

#ifdef DEBUG_CQ
        cerr << "octave " << i << ": resampler latency = " << latencies[i]
             << ", drop " << drops[i] << " (/factor = " << drops[i]/factor
             << "), octaveLatency = " << octaveLatency << " -> "
             << int(round(octaveLatency)) << " (diff * factor = "
             << (octaveLatency - round(octaveLatency)) << " * "
             << factor << " = "
             << (octaveLatency - round(octaveLatency)) * factor << ")" << endl;

        cerr << "double(" << totalLatency << " - " 
             << latencies[i] << " - " << drops[i] << " + " 
             << m_bigBlockSize << ") / " << factor << " = " 
             << octaveLatency << endl;
#endif

        m_buffers.push_back
            (RealSequence(int(octaveLatency + 0.5), 0.0));
    }

    m_fft = new FFTReal(m_p.fftSize);
}

ConstantQ::ComplexBlock
ConstantQ::process(const RealSequence &td)
{
    m_buffers[0].insert(m_buffers[0].end(), td.begin(), td.end());

    for (int i = 1; i < m_octaves; ++i) {
        RealSequence dec = m_decimators[i]->process(td.data(), td.size());
        m_buffers[i].insert(m_buffers[i].end(), dec.begin(), dec.end());
    }

    ComplexBlock out;

    while (true) {

	// We could have quite different remaining sample counts in
	// different octaves, because (apart from the predictable
	// added counts for decimator output on each block) we also
	// have variable additional latency per octave
	bool enough = true;
	for (int i = 0; i < m_octaves; ++i) {
	    int required = m_p.fftSize * pow(2, m_octaves - i - 1);
	    if ((int)m_buffers[i].size() < required) {
		enough = false;
	    }
	}
	if (!enough) break;

        int base = out.size();
        int totalColumns = pow(2, m_octaves - 1) * m_p.atomsPerFrame;
        for (int i = 0; i < totalColumns; ++i) {
            out.push_back(ComplexColumn());
        }

        for (int octave = 0; octave < m_octaves; ++octave) {

            int blocksThisOctave = pow(2, (m_octaves - octave - 1));

            for (int b = 0; b < blocksThisOctave; ++b) {
                ComplexBlock block = processOctaveBlock(octave);
                
                for (int j = 0; j < m_p.atomsPerFrame; ++j) {

                    int target = base +
			    (b * (totalColumns / blocksThisOctave) + 
			     (j * ((totalColumns / blocksThisOctave) /
				   m_p.atomsPerFrame)));

                    while (int(out[target].size()) < 
                           m_p.binsPerOctave * (octave + 1)) {
                        out[target].push_back(Complex());
                    }
                    
                    for (int i = 0; i < m_p.binsPerOctave; ++i) {
                        out[target][m_p.binsPerOctave * octave + i] = 
                            block[j][m_p.binsPerOctave - i - 1];
                    }
                }
            }
        }
    }

    return out;
}

ConstantQ::ComplexBlock
ConstantQ::getRemainingOutput()
{
    // Same as padding added at start, though rounded up
    int pad = ceil(double(m_outputLatency) / m_bigBlockSize) * m_bigBlockSize;
    RealSequence zeros(pad, 0.0);
    return process(zeros);
}

ConstantQ::ComplexBlock
ConstantQ::processOctaveBlock(int octave)
{
    RealSequence ro(m_p.fftSize, 0.0);
    RealSequence io(m_p.fftSize, 0.0);

    m_fft->forward(m_buffers[octave].data(), ro.data(), io.data());

    m_buffers[octave] = RealSequence(m_buffers[octave].begin() + m_p.fftHop,
                                     m_buffers[octave].end());

    ComplexSequence cv(m_p.fftSize);
    for (int i = 0; i < m_p.fftSize; ++i) {
        cv[i] = Complex(ro[i], io[i]);
    }

    ComplexSequence cqrowvec = m_kernel->processForward(cv);

    // Reform into a column matrix
    ComplexBlock cqblock;
    for (int j = 0; j < m_p.atomsPerFrame; ++j) {
        cqblock.push_back(ComplexColumn());
        for (int i = 0; i < m_p.binsPerOctave; ++i) {
            cqblock[j].push_back(cqrowvec[i * m_p.atomsPerFrame + j]);
        }
    }

    return cqblock;
}


