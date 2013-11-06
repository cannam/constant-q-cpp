
#include "ConstantQ.h"

#include "CQKernel.h"

#include "qm-dsp/dsp/rateconversion/Resampler.h"
#include "qm-dsp/maths/MathUtilities.h"
#include "qm-dsp/dsp/transforms/FFT.h"

#include <algorithm>
#include <complex>
#include <iostream>
#include <stdexcept>

using std::vector;
using std::complex;
using std::cerr;
using std::endl;

typedef std::complex<double> C;

ConstantQ::ConstantQ(double sampleRate,
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

ConstantQ::~ConstantQ()
{
    delete m_fft;
    for (int i = 0; i < m_decimators.size(); ++i) {
	delete m_decimators[i];
    }
    delete m_kernel;
}

void
ConstantQ::initialise()
{
    m_octaves = int(ceil(log2(m_maxFrequency / m_minFrequency)));
    double actualMinFreq =
	(m_maxFrequency / pow(2.0, m_octaves)) * pow(2.0, 1.0/m_binsPerOctave);

    cerr << "actual min freq = " << actualMinFreq << endl;

    m_kernel = new CQKernel(m_sampleRate, m_maxFrequency, m_binsPerOctave);
    m_p = m_kernel->getProperties();
    
    // use exact powers of two for resampling rates. They don't have
    // to be related to our actual samplerate, the resampler only
    // cares about the ratio

    int sourceRate = pow(2, m_octaves);
    vector<int> latencies;

    // top octave, no resampling
    latencies.push_back(0);
    m_decimators.push_back(0);

    for (int oct = 1; oct < m_octaves; ++oct) {
	Resampler *r = new Resampler
	    (sourceRate, sourceRate / pow(2, oct), 60, 0.02);
	latencies.push_back(r->getLatency());
	m_decimators.push_back(r);
    }

    //!!! should be multiple of the kernel fft size?
    int maxLatency = *std::max_element(latencies.begin(), latencies.end());
    m_totalLatency = MathUtilities::nextPowerOfTwo(maxLatency);
    cerr << "total latency = " << m_totalLatency << endl;
    for (int i = 0; i < latencies.size(); ++i) {
	m_extraLatencies.push_back(m_totalLatency - latencies[i]);
	cerr << "extra latency " << i << " = " << m_extraLatencies[i] << endl;
	m_buffers.push_back(vector<double>(m_extraLatencies[i], 0.0));
    }

    m_fft = new FFTReal(m_p.fftSize);
    m_bigBlockSize = m_p.fftSize * pow(2, m_octaves) / 2;

    cerr << "m_bigBlockSize = " << m_bigBlockSize << " for " << m_octaves << " octaves" << endl;
}

vector<vector<double> > 
ConstantQ::process(vector<double> td)
{
    m_buffers[0].insert(m_buffers[0].end(), td.begin(), td.end());

    for (int i = 1; i < m_octaves; ++i) {
	vector<double> dec = m_decimators[i]->process(td.data(), td.size());
	m_buffers[i].insert(m_buffers[i].end(), dec.begin(), dec.end());
    }

    vector<vector<double> > out;

    //!!!! need some mechanism for handling remaining samples at the end

    while (m_buffers[0].size() >= m_bigBlockSize) {

	int base = out.size();
	int totalColumns = pow(2, m_octaves - 1) * m_p.atomsPerFrame;
	cerr << "totalColumns = " << totalColumns << endl;
	for (int i = 0; i < totalColumns; ++i) {
	    out.push_back(vector<double>());
	}

	for (int octave = 0; octave < m_octaves; ++octave) {

	    int blocksThisOctave = pow(2, (m_octaves - octave - 1));
//	    cerr << "octave " << octave+1 << " of " << m_octaves << ", n = " << blocksThisOctave << endl;

	    for (int b = 0; b < blocksThisOctave; ++b) {
		vector<vector<double> > block = processOctaveBlock(octave);
		
		for (int j = 0; j < m_p.atomsPerFrame; ++j) {
		    int target = base +
			(b * (totalColumns / blocksThisOctave) + 
			 (j * ((totalColumns / blocksThisOctave) / m_p.atomsPerFrame)));
//		    cerr << "j " << j << " -> target " << target << endl;
		    for (int i = 0; i < m_p.binsPerOctave; ++i) {
			out[target].push_back(block[j][m_p.binsPerOctave-i-1]);
		    }
		}
	    }
	}
    }
    
    return out;
}

vector<vector<double> >
ConstantQ::processOctaveBlock(int octave)
{
    vector<double> ro(m_p.fftSize, 0.0);
    vector<double> io(m_p.fftSize, 0.0);

    m_fft->forward(m_buffers[octave].data(), ro.data(), io.data());

    vector<double> shifted;
    shifted.insert(shifted.end(), 
		   m_buffers[octave].begin() + m_p.fftHop,
		   m_buffers[octave].end());
    m_buffers[octave] = shifted;

    vector<C> cv;
    for (int i = 0; i < m_p.fftSize; ++i) {
	cv.push_back(C(ro[i], io[i]));
    }

    vector<C> cqrowvec = m_kernel->process(cv);

    // Reform into a column matrix 
    vector<vector<double> > cqblock;
    for (int j = 0; j < m_p.atomsPerFrame; ++j) {
	cqblock.push_back(vector<double>());
	for (int i = 0; i < m_p.binsPerOctave; ++i) {
	    cqblock[j].push_back(abs(cqrowvec[i * m_p.atomsPerFrame + j]));
	}
    }

    return cqblock;
}


