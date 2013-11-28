/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#ifndef CQVAMP_H
#define CQVAMP_H

#include <vamp-sdk/Plugin.h>

class ConstantQ;

class CQVamp : public Vamp::Plugin
{
public:
    CQVamp(float inputSampleRate);
    virtual ~CQVamp();

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    InputDomain getInputDomain() const { return TimeDomain; }

    std::string getIdentifier() const;
    std::string getName() const;
    std::string getDescription() const;
    std::string getMaker() const;
    int getPluginVersion() const;
    std::string getCopyright() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(std::string) const;
    void setParameter(std::string, float);

    size_t getPreferredStepSize() const;
    size_t getPreferredBlockSize() const;

    OutputList getOutputDescriptors() const;

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    ConstantQ *m_cq;
    float m_maxFrequency;
    float m_minFrequency;
    int m_bpo;
    int m_stepSize;
    int m_blockSize;

    Vamp::RealTime m_startTime;
    bool m_haveStartTime;
    int m_columnCount;

    std::vector<float> m_prevFeature;
    FeatureSet convertToFeatures(const std::vector<std::vector<double> > &);
};


#endif
