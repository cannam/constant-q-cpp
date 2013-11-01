
#include "CQKernel.h"

#include "qm-dsp/maths/MathUtilities.h"
#include "qm-dsp/dsp/transforms/FFT.h"
#include "qm-dsp/base/Window.h"

#include <cmath>
#include <cassert>
#include <vector>

using std::vector;

CQKernel::CQKernel(double sampleRate, double maxFreq, int binsPerOctave)
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

    int lastCentre = m_p.firstCentre + (m_p.atomsPerFrame - 1) * m_p.atomSpacing;

    m_p.fftHop = (lastCentre + m_p.atomSpacing) - m_p.firstCentre;

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
	
	double fk = m_p.minFrequency * pow(2, ((k-1.0) * bpo));

	vector<double> reals, imags;
	
	for (int i = 0; i < nk; ++i) {
	    double arg = 2.0 * M_PI * fk * i / m_p.sampleRate;
	    reals[i] = win[i] * cos(arg);
	    imags[i] = win[i] * sin(arg);
	}

	int atomOffset = m_p.firstCentre - ceil(nk/2);

	for (int i = 0; i < m_p.atomsPerFrame; ++i) {

	    int shift = atomOffset + (i * m_p.atomSpacing);

	    vector<double> rin(m_p.fftSize, 0.0);
	    vector<double> iin(m_p.fftSize, 0.0);

	    for (int j = 0; j < nk; ++j) {
		rin[j + shift] = reals[j];
		iin[i + shift] = imags[j];
	    }

	    vector<double> rout(m_p.fftSize, 0.0);
	    vector<double> iout(m_p.fftSize, 0.0);
		
	    m_fft->process(false,
			   rin.data(), iin.data(),
			   rout.data(), iout.data());

	    int firstNZ = -1, lastNZ = -1;

	    for (int j = 0; j < m_p.fftSize; ++j) {
		if (sqrt(rout[i] * rout[i] + iout[i] * iout[i]) >= thresh) {
		    lastNZ = j;
		    if (firstNZ < 0) firstNZ = j;
		}
	    }

	    vector<double> rnz, inz;

	    if (firstNZ >= 0) {
		for (int j = firstNZ; j <= lastNZ; ++j) {
		    rnz.push_back(rout[j] / m_p.fftSize);
		    inz.push_back(iout[j] / m_p.fftSize);
		}
		m_kernel.offsets.push_back(firstNZ);
	    } else {
		m_kernel.offsets.push_back(0);
	    }

	    m_kernel.real.push_back(rnz);
	    m_kernel.imag.push_back(inz);
	}
    }

    assert((int)m_kernel.offsets.size() == m_p.binsPerOctave);
    assert((int)m_kernel.real.size() == m_p.binsPerOctave);
    assert((int)m_kernel.imag.size() == m_p.binsPerOctave);

    //!!! and normalise
}


