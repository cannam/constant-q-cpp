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

#include "CQVamp.h"

#include "base/Pitch.h"

#include <algorithm>
#include <cstdio>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

CQVamp::CQVamp(float inputSampleRate) :
    Vamp::Plugin(inputSampleRate),
    m_minMIDIPitch(36),
    m_maxMIDIPitch(84),
    m_tuningFrequency(440),
    m_bpo(24),
    m_interpolation(CQSpectrogram::InterpolateLinear),
    m_cq(0),
    m_maxFrequency(inputSampleRate/2),
    m_minFrequency(46),
    m_haveStartTime(false),
    m_columnCount(0)
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
    desc.identifier = "minpitch";
    desc.name = "Minimum Pitch";
    desc.unit = "MIDI units";
    desc.description = "MIDI pitch corresponding to the lowest frequency to be included in the constant-Q transform";
    desc.minValue = 0;
    desc.maxValue = 127;
    desc.defaultValue = 36;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    desc.identifier = "maxpitch";
    desc.name = "Maximum Pitch";
    desc.unit = "MIDI units";
    desc.description = "MIDI pitch corresponding to the highest frequency to be included in the constant-Q transform";
    desc.minValue = 0;
    desc.maxValue = 127;
    desc.defaultValue = 84;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    desc.identifier = "tuning";
    desc.name = "Tuning Frequency";
    desc.unit = "Hz";
    desc.description = "Frequency of concert A";
    desc.minValue = 360;
    desc.maxValue = 500;
    desc.defaultValue = 440;
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

    desc.identifier = "interpolation";
    desc.name = "Interpolation";
    desc.unit = "";
    desc.description = "Interpolation method used to fill empty cells in lower octaves";
    desc.minValue = 0;
    desc.maxValue = 2;
    desc.defaultValue = 2;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.push_back("None, leave as zero");
    desc.valueNames.push_back("None, repeat prior value");
    desc.valueNames.push_back("Linear interpolation");
    list.push_back(desc);

    return list;
}

float
CQVamp::getParameter(std::string param) const
{
    if (param == "minpitch") {
        return m_minMIDIPitch;
    }
    if (param == "maxpitch") {
        return m_maxMIDIPitch;
    }
    if (param == "tuning") {
        return m_tuningFrequency;
    }
    if (param == "bpo") {
        return m_bpo;
    }
    if (param == "interpolation") {
        return (float)m_interpolation;
    }
    std::cerr << "WARNING: CQVamp::getParameter: unknown parameter \""
              << param << "\"" << std::endl;
    return 0.0;
}

void
CQVamp::setParameter(std::string param, float value)
{
    if (param == "minpitch") {
        m_minMIDIPitch = lrintf(value);
    } else if (param == "maxpitch") {
        m_maxMIDIPitch = lrintf(value);
    } else if (param == "tuning") {
        m_tuningFrequency = value;
    } else if (param == "bpo") {
        m_bpo = lrintf(value);
    } else if (param == "interpolation") {
        m_interpolation = (CQSpectrogram::Interpolation)lrintf(value);
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

    m_minFrequency = Pitch::getFrequencyForPitch
        (m_minMIDIPitch, 0, m_tuningFrequency);
    m_maxFrequency = Pitch::getFrequencyForPitch
        (m_maxMIDIPitch, 0, m_tuningFrequency);

    m_cq = new CQSpectrogram
	(m_inputSampleRate, m_minFrequency, m_maxFrequency, m_bpo,
         m_interpolation);

    return true;
}

void
CQVamp::reset()
{
    if (m_cq) {
	delete m_cq;
	m_cq = new CQSpectrogram
	    (m_inputSampleRate, m_minFrequency, m_maxFrequency, m_bpo,
             m_interpolation);
    }
    m_prevFeature.clear();
    m_haveStartTime = false;
    m_columnCount = 0;
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

    if (m_cq) {
        char name[20];
        for (int i = 0; i < d.binCount; ++i) {
            float freq = m_cq->getBinFrequency(i);
            sprintf(name, "%.1f Hz", freq);
            d.binNames.push_back(name);
        }
    }

    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleType = OutputDescriptor::FixedSampleRate;
    d.sampleRate = m_inputSampleRate / (m_cq ? m_cq->getColumnHop() : 256);
    list.push_back(d);

    return list;
}

CQVamp::FeatureSet
CQVamp::process(const float *const *inputBuffers,
		Vamp::RealTime timestamp)
{
    if (!m_cq) {
	cerr << "ERROR: CQVamp::process: "
	     << "Plugin has not been initialised"
	     << endl;
	return FeatureSet();
    }

    if (!m_haveStartTime) {
        m_startTime = timestamp;
        m_haveStartTime = true;
    }

    vector<double> data;
    for (int i = 0; i < m_blockSize; ++i) data.push_back(inputBuffers[0][i]);
    
    vector<vector<double> > cqout = m_cq->process(data);
    return convertToFeatures(cqout);
}

CQVamp::FeatureSet
CQVamp::getRemainingFeatures()
{
    vector<vector<double> > cqout = m_cq->getRemainingOutput();
    return convertToFeatures(cqout);
}

CQVamp::FeatureSet
CQVamp::convertToFeatures(const vector<vector<double> > &cqout)
{
    FeatureSet returnFeatures;

    int width = cqout.size();
    int height = m_cq->getTotalBins();

    for (int i = 0; i < width; ++i) {

	vector<float> column(height, 0.f);
        int thisHeight = cqout[i].size();
	for (int j = 0; j < thisHeight; ++j) {
	    column[j] = cqout[i][j];
	}

        // put low frequencies at the start
        std::reverse(column.begin(), column.end());

	m_prevFeature = column;

	Feature feature;
	feature.hasTimestamp = true;
        feature.timestamp = m_startTime + Vamp::RealTime::frame2RealTime
            (m_columnCount * m_cq->getColumnHop() - m_cq->getLatency(),
             m_inputSampleRate);
	feature.values = column;
	feature.label = "";

//        cerr << "timestamp = " << feature.timestamp << " (start time = " << m_startTime << ", column count = " << m_columnCount << ", latency = " << m_cq->getLatency() << ", sample rate " << m_inputSampleRate << ")" << endl;

        if (feature.timestamp >= m_startTime) {
            returnFeatures[0].push_back(feature);
        }

        ++m_columnCount;
    }

    return returnFeatures;
}

