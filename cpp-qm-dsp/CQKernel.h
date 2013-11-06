/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#ifndef CQ_KERNEL_H
#define CQ_KERNEL_H

#include <vector>
#include <complex>

class FFT;

class CQKernel
{
public:
    CQKernel(double sampleRate, double maxFreq, int binsPerOctave);
    ~CQKernel();
    
    struct Properties {
        double sampleRate;
        double maxFrequency;
        double minFrequency;
        int binsPerOctave;
        int fftSize;
        int fftHop;
        int atomsPerFrame;
        int atomSpacing;
        int firstCentre;
        double Q;
    };

    Properties getProperties() const { return m_p; }

    std::vector<std::complex<double> > process
        (const std::vector<std::complex<double> > &);

private:
    Properties m_p;
    FFT *m_fft;

    struct KernelMatrix {
        std::vector<int> origin;
        std::vector<std::vector<std::complex<double> > > data;
    };
    KernelMatrix m_kernel;

    void generateKernel();
    void finaliseKernel();
};

#endif
