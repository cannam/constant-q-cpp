/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "CQVamp.h"

#include "../cpp-qm-dsp/ConstantQ.h"

using std::string;
using std::vector;
using std::cerr;
using std::endl;

CQVamp::CQVamp(float inputSampleRate) :
    Vamp::Plugin(inputSampleRate),
    m_cq(0),
    m_maxFrequency(inputSampleRate/2),
    m_minFrequency(46),
    m_bpo(24)
{
}

CQVamp::~CQVamp()
{
    delete m_cq;
}

string
CQVamp::getIdentifier() const
{
    return "cqvamp";
}

string
CQVamp::getName() const
{
    return "Constant-Q Spectrogram";
}

string
CQVamp::getDescription() const
{
    return "Extract a spectrogram with constant ratio of centre frequency to resolution from the input audio";
}

string
CQVamp::getMaker() const
{
    return "Queen Mary, University of London";
}

int
CQVamp::getPluginVersion() const
{
    return 1;
}

string
CQVamp::getCopyright() const
{
    return "Plugin by Chris Cannam. Method by Christian Schörkhuber and Anssi Klapuri. Copyright (c) 2013 QMUL";
}

CQVamp::ParameterList
CQVamp::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor desc;
    desc.identifier = "minfreq";
    desc.name = "Minimum Frequency";
    desc.unit = "Hz";
    desc.description = "Hint for the lowest frequency to be included in the constant-Q transform. The actual frequency range will be an integral number of octaves ending at the highest frequency specified";
    desc.minValue = 10;
    desc.maxValue = m_inputSampleRate/2;
    desc.defaultValue = 46;
    desc.isQuantized = false;
    list.push_back(desc);

    desc.identifier = "maxfreq";
    desc.name = "Maximum Frequency";
    desc.unit = "Hz";
    desc.description = "Highest frequency to be included in the constant-Q transform";
    desc.minValue = 10;
    desc.maxValue = m_inputSampleRate/2;
    desc.defaultValue = m_inputSampleRate/2;
    desc.isQuantized = false;
    list.push_back(desc);
    
    desc.identifier = "bpo";
    desc.name = "Bins per Octave";
    desc.unit = "bins";
    desc.description = "Number of constant-Q transform bins per octave";
    desc.minValue = 2;
    desc.maxValue = 480;
    desc.defaultValue = 24;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    return list;
}

float
CQVamp::getParameter(std::string param) const
{
    if (param == "minfreq") {
        return m_minFrequency;
    }
    if (param == "maxfreq") {
        return m_maxFrequency;
    }
    if (param == "bpo") {
        return m_bpo;
    }
    std::cerr << "WARNING: CQVamp::getParameter: unknown parameter \""
              << param << "\"" << std::endl;
    return 0.0;
}

void
CQVamp::setParameter(std::string param, float value)
{
    if (param == "minfreq") {
        m_minFrequency = value;
    } else if (param == "maxfreq") {
        m_maxFrequency = value;
    } else if (param == "bpo") {
        m_bpo = lrintf(value);
    } else {
        std::cerr << "WARNING: CQVamp::setParameter: unknown parameter \""
                  << param << "\"" << std::endl;
    }
}

bool
CQVamp::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (m_cq) {
	delete m_cq;
	m_cq = 0;
    }

    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

    m_stepSize = stepSize;
    m_blockSize = blockSize;

    m_cq = new ConstantQ
	(m_inputSampleRate, m_minFrequency, m_maxFrequency, m_bpo);

    return true;
}

void
CQVamp::reset()
{
    if (m_cq) {
	delete m_cq;
	m_cq = new ConstantQ
	    (m_inputSampleRate, m_minFrequency, m_maxFrequency, m_bpo);
    }
    m_prevFeature.clear();
}

size_t
CQVamp::getPreferredStepSize() const
{
    return 0;
}

size_t
CQVamp::getPreferredBlockSize() const
{
    return 0;
}

CQVamp::OutputList
CQVamp::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor d;
    d.identifier = "constantq";
    d.name = "Constant-Q Spectrogram";
    d.unit = "";
    d.description = "Output of constant-Q transform, as a single vector per process block";
    d.hasFixedBinCount = true;
    d.binCount = (m_cq ? m_cq->getTotalBins() : (9 * 24));
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = m_inputSampleRate / (m_cq ? m_cq->getColumnHop() : 256);
    list.push_back(d);

    return list;
}

CQVamp::FeatureSet
CQVamp::process(const float *const *inputBuffers,
		Vamp::RealTime /* timestamp */)
{
    if (!m_cq) {
	cerr << "ERROR: CQVamp::process: "
	     << "Plugin has not been initialised"
	     << endl;
	return FeatureSet();
    }

    vector<double> data;
    for (int i = 0; i < m_blockSize; ++i) data.push_back(inputBuffers[0][i]);
    
    vector<vector<double> > cqout = m_cq->process(data);
    return convertToFeatures(cqout);
}

CQVamp::FeatureSet
CQVamp::getRemainingFeatures()
{
    vector<vector<double> > cqout = m_cq->getRemainingBlocks();
    return convertToFeatures(cqout);
}

CQVamp::FeatureSet
CQVamp::convertToFeatures(const vector<vector<double> > &cqout)
{
    
    FeatureSet returnFeatures;

    for (int i = 0; i < (int)cqout.size(); ++i) {

	vector<float> column(m_cq->getTotalBins(), 0.f);

	for (int j = 0; j < (int)cqout[i].size(); ++j) {
	    column[j] = cqout[i][j];
	}
	for (int j = cqout[i].size(); j < m_cq->getTotalBins(); ++j) {
	    if (j < (int)m_prevFeature.size()) {
		column[j] = m_prevFeature[j];
	    }
	}

	m_prevFeature = column;

	Feature feature;
	feature.hasTimestamp = false;
	feature.values = column;
	feature.label = "";
	returnFeatures[0].push_back(feature);
    }

    return returnFeatures;
}

