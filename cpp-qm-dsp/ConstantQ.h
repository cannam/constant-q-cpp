/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#ifndef CONSTANTQ_H
#define CONSTANTQ_H

#include "CQKernel.h"

#include <vector>

class Resampler;
class FFTReal;

class ConstantQ
{
public:
    ConstantQ(double sampleRate, 
	      double minFreq, double maxFreq, 
	      int binsPerOctave);
    ~ConstantQ();

    double getSampleRate() const { return m_sampleRate; }
    double getMaxFrequency() const { return m_p.maxFrequency; }
    double getMinFrequency() const { return m_p.minFrequency; } // actual min, not that provided to ctor
    int getBinsPerOctave() const { return m_binsPerOctave; }
    int getOctaves() const { return m_octaves; }
    int getTotalBins() const { return m_octaves * m_binsPerOctave; }
    int getColumnHop() const { return m_p.fftHop / m_p.atomsPerFrame; }

    std::vector<std::vector<double> > process(const std::vector<double> &);
    std::vector<std::vector<double> > getRemainingBlocks();

private:
    double m_sampleRate;
    double m_maxFrequency;
    double m_minFrequency;
    int m_binsPerOctave;
    int m_octaves;

    CQKernel *m_kernel;
    CQKernel::Properties m_p;
    int m_bigBlockSize;

    std::vector<Resampler *> m_decimators;
    std::vector<std::vector<double> > m_buffers;

    int m_totalLatency;

    FFTReal *m_fft;

    void initialise();
    std::vector<std::vector<double> > processOctaveBlock(int octave);
};

#endif

