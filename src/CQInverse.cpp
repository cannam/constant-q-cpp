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

CQInverse::CQInverse(CQParameters params) :
    m_inparams(params),
    m_sampleRate(params.sampleRate),
    m_maxFrequency(params.maxFrequency),
    m_minFrequency(params.minFrequency),
    m_binsPerOctave(params.binsPerOctave),
    m_fft(0)
{
    if (m_minFrequency <= 0.0 || m_maxFrequency <= 0.0) {
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
CQInverse::getBinFrequency(double bin) const
{
    // our bins are returned in high->low order
    bin = (getBinsPerOctave() * getOctaves()) - bin - 1;
    return getMinFrequency() * pow(2, (bin / getBinsPerOctave()));
}

void
CQInverse::initialise()
{
    m_octaves = int(ceil(log(m_maxFrequency / m_minFrequency) / log(2)));

    if (m_octaves < 1) {
        m_kernel = 0; // incidentally causing isValid() to return false
        return;
    }

    m_kernel = new CQKernel(m_inparams);
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
            (sourceRate / factor, sourceRate, 50, 0.05);

#ifdef DEBUG_CQ
        cerr << "inverse: octave " << i << ": resample from " << sourceRate/factor << " to " << sourceRate << endl;
#endif

	// See ConstantQ.cpp for discussion on latency -- output
	// latency here is at target rate which, this way around, is
	// what we want

        latencies.push_back(r->getLatency());
        m_upsamplers.push_back(r);
    }

    // additionally we will have fftHop latency at individual octave
    // rate (before upsampling) for the overlap-add in each octave
    for (int i = 0; i < m_octaves; ++i) {
        latencies[i] += m_p.fftHop * pow(2, i);
    }

    // Now reverse the drop adjustment made in ConstantQ to align the
    // atom centres across different octaves (but this time at output
    // sample rate)

    int emptyHops = m_p.firstCentre / m_p.atomSpacing;

    vector<int> pushes;
    for (int i = 0; i < m_octaves; ++i) {
	int factor = pow(2, i);
	int pushHops = emptyHops * pow(2, m_octaves - i - 1) - emptyHops;
	int push = ((pushHops * m_p.fftHop) * factor) / m_p.atomsPerFrame;
	pushes.push_back(push);
    }

    int maxLatLessPush = 0;
    for (int i = 0; i < m_octaves; ++i) {
	int latLessPush = latencies[i] - pushes[i];
	if (latLessPush > maxLatLessPush) maxLatLessPush = latLessPush;
    }

    int totalLatency = maxLatLessPush + 10;
    if (totalLatency < 0) totalLatency = 0;

    m_outputLatency = totalLatency + m_p.firstCentre * pow(2, m_octaves-1);

#ifdef DEBUG_CQ
    cerr << "totalLatency = " << totalLatency << ", m_outputLatency = " << m_outputLatency << endl;
#endif

    for (int i = 0; i < m_octaves; ++i) {

	// Calculate the difference between the total latency applied
	// across all octaves, and the existing latency due to the
	// upsampler for this octave.

        int latencyPadding = totalLatency - latencies[i] + pushes[i];

#ifdef DEBUG_CQ
        cerr << "octave " << i << ": push " << pushes[i] << ", resampler latency inc overlap space " << latencies[i] << ", latencyPadding = " << latencyPadding << " (/factor = " << latencyPadding / pow(2, i) << ")" << endl;
#endif

        m_buffers.push_back(RealSequence(latencyPadding, 0.0));
    }

    for (int i = 0; i < m_octaves; ++i) {
        // Fixed-size buffer for IFFT overlap-add
        m_olaBufs.push_back(RealSequence(m_p.fftSize, 0.0));
    }

    m_fft = new FFTReal(m_p.fftSize);
}

CQInverse::RealSequence
CQInverse::process(const ComplexBlock &block)
{
    // The input data is of the form produced by ConstantQ::process --
    // an unknown number N of columns of varying height. We assert
    // that N is a multiple of atomsPerFrame * 2^(octaves-1), as must
    // be the case for data that came directly from our ConstantQ
    // implementation.

    int widthProvided = block.size();

    if (widthProvided == 0) {
        return drawFromBuffers();
    }

    int blockWidth = m_p.atomsPerFrame * int(pow(2, m_octaves - 1));

    if (widthProvided % blockWidth != 0) {
        cerr << "ERROR: CQInverse::process: Input block size ("
             << widthProvided
             << ") must be a multiple of processing block width "
             << "(atoms-per-frame * 2^(octaves-1) = "
             << m_p.atomsPerFrame << " * 2^(" << m_octaves << "-1) = "
             << blockWidth << ")" << endl;
        throw std::invalid_argument
            ("Input block size must be a multiple of processing block width");
    }

    // Procedure:
    // 
    // 1. Slice the list of columns into a set of lists of columns,
    // one per octave, each of width N / (2^octave-1) and height
    // binsPerOctave, containing the values present in that octave
    //
    // 2. Group each octave list by atomsPerFrame columns at a time,
    // and stack these so as to achieve a list, for each octave, of
    // taller columns of height binsPerOctave * atomsPerFrame
    //
    // 3. For each taller column, take the product with the inverse CQ
    // kernel (which is the conjugate of the forward kernel) and
    // perform an inverse FFT
    //
    // 4. Overlap-add each octave's resynthesised blocks (unwindowed)
    //
    // 5. Resample each octave's overlap-add stream to the original
    // rate
    //
    // 6. Sum the resampled streams and return
    
    for (int i = 0; i < m_octaves; ++i) {
        
        // Step 1

        ComplexBlock oct;

        for (int j = 0; j < widthProvided; ++j) {
            int h = block[j].size();
            if (h < m_binsPerOctave * (i+1)) {
                continue;
            }
            ComplexColumn col(block[j].begin() + m_binsPerOctave * i,
                              block[j].begin() + m_binsPerOctave * (i+1));
            oct.push_back(col);
        }

        // Steps 2, 3, 4, 5
        processOctave(i, oct);
    }
    
    // Step 6
    return drawFromBuffers();
}

CQInverse::RealSequence
CQInverse::drawFromBuffers()
{
    // 6. Sum the resampled streams and return

    int available = 0;

    for (int i = 0; i < m_octaves; ++i) {
        if (i == 0 || int(m_buffers[i].size()) < available) {
            available = m_buffers[i].size();
        }
    }

    RealSequence result(available, 0);

    if (available == 0) {
        return result;
    }

    for (int i = 0; i < m_octaves; ++i) {
        for (int j = 0; j < available; ++j) {
            result[j] += m_buffers[i][j];
        }
        m_buffers[i] = RealSequence(m_buffers[i].begin() + available,
                                    m_buffers[i].end());
    }

    return result;
}

CQInverse::RealSequence
CQInverse::getRemainingOutput()
{
    for (int j = 0; j < m_octaves; ++j) {
        int factor = pow(2, j);
        int latency = (j > 0 ? m_upsamplers[j]->getLatency() : 0) / factor;
        for (int i = 0; i < (latency + m_p.fftSize) / m_p.fftHop; ++i) {
            overlapAddAndResample(j, RealSequence(m_olaBufs[j].size(), 0));
        }
    }

    return drawFromBuffers();
}

void
CQInverse::processOctave(int octave, const ComplexBlock &columns)
{
    // 2. Group each octave list by atomsPerFrame columns at a time,
    // and stack these so as to achieve a list, for each octave, of
    // taller columns of height binsPerOctave * atomsPerFrame

    int ncols = columns.size();

    if (ncols % m_p.atomsPerFrame != 0) {
        cerr << "ERROR: CQInverse::process: Number of columns ("
             << ncols
             << ") in octave " << octave
             << " must be a multiple of atoms-per-frame ("
             << m_p.atomsPerFrame << ")" << endl;
        throw std::invalid_argument
            ("Columns in octave must be a multiple of atoms per frame");
    }

    for (int i = 0; i < ncols; i += m_p.atomsPerFrame) {

        ComplexColumn tallcol;
        for (int b = 0; b < m_binsPerOctave; ++b) {
            for (int a = 0; a < m_p.atomsPerFrame; ++a) {
                tallcol.push_back(columns[i + a][m_binsPerOctave - b - 1]);
            }
        }
        
        processOctaveColumn(octave, tallcol);
    }
}

void
CQInverse::processOctaveColumn(int octave, const ComplexColumn &column)
{
    // 3. For each taller column, take the product with the inverse CQ
    // kernel (which is the conjugate of the forward kernel) and
    // perform an inverse FFT

    if ((int)column.size() != m_p.atomsPerFrame * m_binsPerOctave) {
        cerr << "ERROR: CQInverse::processOctaveColumn: Height of column ("
             << column.size() << ") in octave " << octave
             << " must be atoms-per-frame * bins-per-octave ("
             << m_p.atomsPerFrame << " * " << m_binsPerOctave << " = "
             << m_p.atomsPerFrame * m_binsPerOctave << ")" << endl;
        throw std::invalid_argument
            ("Column height must match atoms-per-frame * bins-per-octave");
    }

    ComplexSequence transformed = m_kernel->processInverse(column);

    int halfLen = m_p.fftSize/2 + 1;

    RealSequence ri(halfLen, 0);
    RealSequence ii(halfLen, 0);

    for (int i = 0; i < halfLen; ++i) {
        ri[i] = transformed[i].real();
        ii[i] = transformed[i].imag();
    }

    RealSequence timeDomain(m_p.fftSize, 0);

    m_fft->inverse(ri.data(), ii.data(), timeDomain.data());

    overlapAddAndResample(octave, timeDomain);
}

void
CQInverse::overlapAddAndResample(int octave, const RealSequence &seq)
{
    // 4. Overlap-add each octave's resynthesised blocks (unwindowed)
    //
    // and
    //
    // 5. Resample each octave's overlap-add stream to the original
    // rate

    if (seq.size() != m_olaBufs[octave].size()) {
        cerr << "ERROR: CQInverse::overlapAdd: input sequence length ("
             << seq.size() << ") is expected to match OLA buffer size ("
             << m_olaBufs[octave].size() << ")" << endl;
        throw std::invalid_argument
            ("Input sequence length should match OLA buffer size");
    }

    RealSequence toResample(m_olaBufs[octave].begin(),
                            m_olaBufs[octave].begin() + m_p.fftHop);

    RealSequence resampled = 
        octave > 0 ?
        m_upsamplers[octave]->process(toResample.data(), toResample.size()) :
        toResample;

    m_buffers[octave].insert(m_buffers[octave].end(),
                             resampled.begin(),
                             resampled.end());
    
    m_olaBufs[octave] = RealSequence(m_olaBufs[octave].begin() + m_p.fftHop,
                                     m_olaBufs[octave].end());
    
    RealSequence pad(m_p.fftHop, 0);

    m_olaBufs[octave].insert(m_olaBufs[octave].end(),
                             pad.begin(),
                             pad.end());

    for (int i = 0; i < m_p.fftSize; ++i) {
        m_olaBufs[octave][i] += seq[i];
    }
}

