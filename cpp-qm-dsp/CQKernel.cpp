/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "CQKernel.h"

#include "qm-dsp/maths/MathUtilities.h"
#include "qm-dsp/dsp/transforms/FFT.h"
#include "qm-dsp/base/Window.h"

#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <algorithm>

using std::vector;
using std::complex;
using std::cerr;
using std::endl;

typedef std::complex<double> C;

CQKernel::CQKernel(double sampleRate, double maxFreq, int binsPerOctave) :
    m_fft(0)
{
    m_p.sampleRate = sampleRate;
    m_p.maxFrequency = maxFreq;
    m_p.binsPerOctave = binsPerOctave;
    generateKernel();
}

CQKernel::~CQKernel()
{
    delete m_fft;
}

void
CQKernel::generateKernel()
{
    double q = 1;
    double atomHopFactor = 0.25;
    double thresh = 0.0005;

    double bpo = m_p.binsPerOctave;

    m_p.minFrequency = (m_p.maxFrequency / 2) * pow(2, 1.0/bpo);
    m_p.Q = q / (pow(2, 1.0/bpo) - 1.0);

    double maxNK = round(m_p.Q * m_p.sampleRate / m_p.minFrequency);
    double minNK = round
        (m_p.Q * m_p.sampleRate /
         (m_p.minFrequency * pow(2, (bpo - 1.0) / bpo)));

    m_p.atomSpacing = round(minNK * atomHopFactor);
    m_p.firstCentre = m_p.atomSpacing * ceil(ceil(maxNK / 2.0) / m_p.atomSpacing);
    m_p.fftSize = MathUtilities::nextPowerOfTwo
        (m_p.firstCentre + ceil(maxNK / 2.0));

    m_p.atomsPerFrame = floor
        (1.0 + (m_p.fftSize - ceil(maxNK / 2.0) - m_p.firstCentre) / m_p.atomSpacing);

    cerr << "atomsPerFrame = " << m_p.atomsPerFrame << " (atomHopFactor = " << atomHopFactor << ")" << endl;

    int lastCentre = m_p.firstCentre + (m_p.atomsPerFrame - 1) * m_p.atomSpacing;

    m_p.fftHop = (lastCentre + m_p.atomSpacing) - m_p.firstCentre;

    cerr << "fftHop = " << m_p.fftHop << endl;

    m_fft = new FFT(m_p.fftSize);

    for (int k = 1; k <= m_p.binsPerOctave; ++k) {
        
        int nk = round(m_p.Q * m_p.sampleRate /
                       (m_p.minFrequency * pow(2, ((k-1.0) / bpo))));

        // The MATLAB version uses a symmetric window, but our windows
        // are periodic. A symmetric window of size N is a periodic
        // one of size N-1 with the first element stuck on the end
        Window<double> w(BlackmanHarrisWindow, nk-1);
        vector<double> win = w.getWindowData();
        win.push_back(win[0]);

        for (int i = 0; i < (int)win.size(); ++i) {
            win[i] = sqrt(win[i]) / nk;
        }

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

    cerr << "size = " << m_kernel.data.size() << "*" << m_kernel.data[0].size() << " (fft size = " << m_p.fftSize << ")" << endl;

    assert((int)m_kernel.data.size() == m_p.binsPerOctave * m_p.atomsPerFrame);
    assert((int)m_kernel.data[0].size() == m_p.fftSize);

    cerr << "density = " << double(nnz) / double(m_p.binsPerOctave * m_p.atomsPerFrame * m_p.fftSize) << " (" << nnz << " of " << m_p.binsPerOctave * m_p.atomsPerFrame * m_p.fftSize << ")" << endl;

    finaliseKernel();
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
    double q = 1; //!!! duplicated from constructor
    for (int i = round(1.0/q); i < ncols - round(1.0/q) - 2; ++i) {
        wK.push_back(abs(square[i][i]));
    }

    double weight = double(m_p.fftHop) / m_p.fftSize;
    weight /= MathUtilities::mean(wK.data(), wK.size());
    weight = sqrt(weight);
    
    cerr << "weight = " << weight << endl;

    // apply normalisation weight, make sparse, and store conjugates
    // (our multiplication order means we will effectively be using
    // the adjoint or conjugate transpose of the kernel matrix)

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
CQKernel::process(const vector<C> &cv)
{
    // matrix multiply m_kernel.data by in, converting in to complex
    // as we go

    int nrows = m_p.binsPerOctave * m_p.atomsPerFrame;

    vector<C> rv(nrows, C(0, 0));

    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < (int)m_kernel.data[i].size(); ++j) {
            rv[i] += cv[j + m_kernel.origin[i]] * m_kernel.data[i][j];
        }
    }

    return rv;
}


