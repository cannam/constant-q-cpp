/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#ifndef CONSTANTQ_H
#define CONSTANTQ_H

#include "CQKernel.h"

#include <vector>

class Resampler;
class FFT;

class ConstantQ
{
public:
    ConstantQ(double sampleRate, 
	      double minFreq, double maxFreq, 
	      int binsPerOctave);
    ~ConstantQ();

    std::vector<std::vector<double> > process(std::vector<double>);

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
    std::vector<std::vector<double> > m_octaveBuffers;

    int m_totalLatency;
    std::vector<int> m_extraLatencies; // per resampler, to make up to total

    FFT *m_fft;

    void initialise();
};

#endif

