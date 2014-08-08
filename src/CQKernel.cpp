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

#include "CQKernel.h"

#include "dsp/MathUtilities.h"
#include "dsp/FFT.h"
#include "dsp/Window.h"

#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <algorithm>

#include <cmath>

using std::vector;
using std::complex;
using std::cerr;
using std::endl;

typedef std::complex<double> C;

//#define DEBUG_CQ_KERNEL 1

CQKernel::CQKernel(CQParameters params) :
    m_inparams(params),
    m_valid(false),
    m_fft(0)
{
    m_p.sampleRate = params.sampleRate;
    m_p.maxFrequency = params.maxFrequency;
    m_p.binsPerOctave = params.binsPerOctave;
    m_valid = generateKernel();
}

CQKernel::~CQKernel()
{
    delete m_fft;
}

vector<double>
CQKernel::makeWindow(int len) const
{
    // The MATLAB version uses a symmetric window, but our windows
    // are periodic. A symmetric window of size N is a periodic
    // one of size N-1 with the first element stuck on the end.

    WindowType wt(BlackmanHarrisWindow);

    switch (m_inparams.window) {
    case CQParameters::SqrtBlackmanHarris:
    case CQParameters::BlackmanHarris:
        wt = BlackmanHarrisWindow;
        break;
    case CQParameters::SqrtBlackman:
    case CQParameters::Blackman:
        wt = BlackmanWindow;
        break;
    case CQParameters::SqrtHann:
    case CQParameters::Hann:
        wt = HanningWindow;
        break;
    }

    Window<double> w(wt, len-1);
    vector<double> win = w.getWindowData();
    win.push_back(win[0]);

    switch (m_inparams.window) {
    case CQParameters::SqrtBlackmanHarris:
    case CQParameters::SqrtBlackman:
    case CQParameters::SqrtHann:
        for (int i = 0; i < (int)win.size(); ++i) {
            win[i] = sqrt(win[i]) / len;
        }
        break;
    case CQParameters::BlackmanHarris:
    case CQParameters::Blackman:
    case CQParameters::Hann:
        for (int i = 0; i < (int)win.size(); ++i) {
            win[i] = win[i] / len;
        }
        break;
    }

    return win;
}

bool
CQKernel::generateKernel()
{
    double q = m_inparams.q;
    double atomHopFactor = m_inparams.atomHopFactor;
    double thresh = m_inparams.threshold;

    double bpo = m_p.binsPerOctave;

    m_p.minFrequency = (m_p.maxFrequency / 2) * pow(2, 1.0/bpo);
    m_p.Q = q / (pow(2, 1.0/bpo) - 1.0);

    double maxNK = int(m_p.Q * m_p.sampleRate / m_p.minFrequency + 0.5);
    double minNK = int
        (m_p.Q * m_p.sampleRate /
         (m_p.minFrequency * pow(2, (bpo - 1.0) / bpo)) + 0.5);

    if (minNK == 0 || maxNK == 0) {
        // most likely pathological parameters of some sort
        cerr << "WARNING: CQKernel::generateKernel: minNK or maxNK is zero (minNK == " << minNK << ", maxNK == " << maxNK << "), not generating a kernel" << endl;
        m_p.atomSpacing = 0;
        m_p.firstCentre = 0;
        m_p.fftSize = 0;
        m_p.atomsPerFrame = 0;
        m_p.lastCentre = 0;
        m_p.fftHop = 0;
        return false;
    }

    m_p.atomSpacing = int(minNK * atomHopFactor + 0.5);
    m_p.firstCentre = m_p.atomSpacing * ceil(ceil(maxNK / 2.0) / m_p.atomSpacing);
    m_p.fftSize = MathUtilities::nextPowerOfTwo
        (m_p.firstCentre + ceil(maxNK / 2.0));

    m_p.atomsPerFrame = floor
        (1.0 + (m_p.fftSize - ceil(maxNK / 2.0) - m_p.firstCentre) / m_p.atomSpacing);

#ifdef DEBUG_CQ_KERNEL
    cerr << "atomsPerFrame = " << m_p.atomsPerFrame << " (q = " << q << ", Q = " << m_p.Q << ", atomHopFactor = " << atomHopFactor << ", atomSpacing = " << m_p.atomSpacing << ", fftSize = " << m_p.fftSize << ", maxNK = " << maxNK << ", firstCentre = " << m_p.firstCentre << ")" << endl;
#endif

    m_p.lastCentre = m_p.firstCentre + (m_p.atomsPerFrame - 1) * m_p.atomSpacing;

    m_p.fftHop = (m_p.lastCentre + m_p.atomSpacing) - m_p.firstCentre;

#ifdef DEBUG_CQ_KERNEL
    cerr << "fftHop = " << m_p.fftHop << endl;
#endif

    m_fft = new FFT(m_p.fftSize);

    for (int k = 1; k <= m_p.binsPerOctave; ++k) {
        
        int nk = int(m_p.Q * m_p.sampleRate /
                     (m_p.minFrequency * pow(2, ((k-1.0) / bpo))) + 0.5);

        vector<double> win = makeWindow(nk);

        double fk = m_p.minFrequency * pow(2, ((k-1.0) / bpo));

        vector<double> reals, imags;
        
        for (int i = 0; i < nk; ++i) {
            double arg = (2.0 * M_PI * fk * i) / m_p.sampleRate;
            reals.push_back(win[i] * cos(arg));
            imags.push_back(win[i] * sin(arg));
        }

        int atomOffset = m_p.firstCentre - int(ceil(nk/2.0));

        for (int i = 0; i < m_p.atomsPerFrame; ++i) {

            int shift = atomOffset + (i * m_p.atomSpacing);

            vector<double> rin(m_p.fftSize, 0.0);
            vector<double> iin(m_p.fftSize, 0.0);

            for (int j = 0; j < nk; ++j) {
                rin[j + shift] = reals[j];
                iin[j + shift] = imags[j];
            }

            vector<double> rout(m_p.fftSize, 0.0);
            vector<double> iout(m_p.fftSize, 0.0);

            m_fft->process(false,
                           rin.data(), iin.data(),
                           rout.data(), iout.data());

            // Keep this dense for the moment (until after
            // normalisation calculations)

            vector<C> row;

            for (int j = 0; j < m_p.fftSize; ++j) {
                if (sqrt(rout[j] * rout[j] + iout[j] * iout[j]) < thresh) {
                    row.push_back(C(0, 0));
                } else {
                    row.push_back(C(rout[j] / m_p.fftSize,
                                    iout[j] / m_p.fftSize));
                }
            }

            m_kernel.origin.push_back(0);
            m_kernel.data.push_back(row);
        }
    }

    assert((int)m_kernel.data.size() == m_p.binsPerOctave * m_p.atomsPerFrame);

    // print density as diagnostic

    int nnz = 0;
    for (int i = 0; i < (int)m_kernel.data.size(); ++i) {
        for (int j = 0; j < (int)m_kernel.data[i].size(); ++j) {
            if (m_kernel.data[i][j] != C(0, 0)) {
                ++nnz;
            }
        }
    }

#ifdef DEBUG_CQ_KERNEL
    cerr << "size = " << m_kernel.data.size() << "*" << m_kernel.data[0].size() << " (fft size = " << m_p.fftSize << ")" << endl;
#endif

    assert((int)m_kernel.data.size() == m_p.binsPerOctave * m_p.atomsPerFrame);
    assert((int)m_kernel.data[0].size() == m_p.fftSize);

#ifdef DEBUG_CQ_KERNEL
    cerr << "density = " << double(nnz) / double(m_p.binsPerOctave * m_p.atomsPerFrame * m_p.fftSize) << " (" << nnz << " of " << m_p.binsPerOctave * m_p.atomsPerFrame * m_p.fftSize << ")" << endl;
#endif

    finaliseKernel();
    return true;
}

static bool ccomparator(C &c1, C &c2)
{
    return abs(c1) < abs(c2);
}

static int maxidx(vector<C> &v)
{
    return std::max_element(v.begin(), v.end(), ccomparator) - v.begin();
}

void
CQKernel::finaliseKernel()
{
    // calculate weight for normalisation

    int wx1 = maxidx(m_kernel.data[0]);
    int wx2 = maxidx(m_kernel.data[m_kernel.data.size()-1]);

    vector<vector<C> > subset(m_kernel.data.size());
    for (int j = wx1; j <= wx2; ++j) {
        for (int i = 0; i < (int)m_kernel.data.size(); ++i) {
            subset[i].push_back(m_kernel.data[i][j]);
        }
    }

    int nrows = subset.size();
    int ncols = subset[0].size();
    vector<vector<C> > square(ncols); // conjugate transpose of subset * subset

    for (int i = 0; i < nrows; ++i) {
        assert((int)subset[i].size() == ncols);
    }

    for (int j = 0; j < ncols; ++j) {
        for (int i = 0; i < ncols; ++i) {
            C v(0, 0);
            for (int k = 0; k < nrows; ++k) {
                v += subset[k][i] * conj(subset[k][j]);
            }
            square[i].push_back(v);
        }
    }

    vector<double> wK;
    double q = m_inparams.q;
    for (int i = int(1.0/q + 0.5); i < ncols - int(1.0/q + 0.5) - 2; ++i) {
        wK.push_back(abs(square[i][i]));
    }

    double weight = double(m_p.fftHop) / m_p.fftSize;
    if (!wK.empty()) {
        weight /= MathUtilities::mean(wK.data(), wK.size());
    }
    weight = sqrt(weight);

#ifdef DEBUG_CQ_KERNEL    
    cerr << "weight = " << weight << " (from " << wK.size() << " elements in wK, ncols = " << ncols << ", q = " << q << ")" << endl;
#endif

    // apply normalisation weight, make sparse, and store conjugate
    // (we use the adjoint or conjugate transpose of the kernel matrix
    // for the forward transform, the plain kernel for the inverse
    // which we expect to be less common)

    KernelMatrix sk;

    for (int i = 0; i < (int)m_kernel.data.size(); ++i) {

        sk.origin.push_back(0);
        sk.data.push_back(vector<C>());

        int lastNZ = 0;
        for (int j = (int)m_kernel.data[i].size()-1; j >= 0; --j) {
            if (abs(m_kernel.data[i][j]) != 0.0) {
                lastNZ = j;
                break;
            }
        }

        bool haveNZ = false;
        for (int j = 0; j <= lastNZ; ++j) {
            if (haveNZ || abs(m_kernel.data[i][j]) != 0.0) {
                if (!haveNZ) sk.origin[i] = j;
                haveNZ = true;
                sk.data[i].push_back(conj(m_kernel.data[i][j]) * weight);
            }
        }
    }

    m_kernel = sk;
}

vector<C>
CQKernel::processForward(const vector<C> &cv)
{
    // straightforward matrix multiply (taking into account m_kernel's
    // slightly-sparse representation)

    if (m_kernel.data.empty()) return vector<C>();

    int nrows = m_p.binsPerOctave * m_p.atomsPerFrame;

    vector<C> rv(nrows, C());

    for (int i = 0; i < nrows; ++i) {
        int len = m_kernel.data[i].size();
        for (int j = 0; j < len; ++j) {
            rv[i] += cv[j + m_kernel.origin[i]] * m_kernel.data[i][j];
        }
    }

    return rv;
}

vector<C>
CQKernel::processInverse(const vector<C> &cv)
{
    // matrix multiply by conjugate transpose of m_kernel. This is
    // actually the original kernel as calculated, we just stored the
    // conjugate-transpose of the kernel because we expect to be doing
    // more forward transforms than inverse ones.

    if (m_kernel.data.empty()) return vector<C>();

    int ncols = m_p.binsPerOctave * m_p.atomsPerFrame;
    int nrows = m_p.fftSize;

    vector<C> rv(nrows, C());

    for (int j = 0; j < ncols; ++j) {
        int i0 = m_kernel.origin[j];
        int i1 = i0 + m_kernel.data[j].size();
        for (int i = i0; i < i1; ++i) {
            rv[i] += cv[j] * conj(m_kernel.data[j][i - i0]);
        }
    }

    return rv;
}


