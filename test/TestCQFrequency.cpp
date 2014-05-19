/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "cq/CQSpectrogram.h"

#include "dsp/Window.h"

#include <cmath>
#include <vector>

using std::vector;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(TestCQFrequency)

// The principle here is to feed a single windowed sinusoid into a
// small CQ transform and check that the output has its peak bin at
// the correct frequency. We can repeat for different frequencies both
// inside and outside the frequency range supported by the CQ. We
// should also repeat for CQSpectrogram outputs as well as the raw CQ.

// Set up fs/2 = 50, frequency range 10 -> 40 i.e. 2 octaves, fixed
// duration of 2 seconds
static const double sampleRate = 100;
static const double cqmin = 10;
static const double cqmax = 40;
static const double bpo = 4;
static const int duration = sampleRate * 2;

// Fairly arbitrary max value for CQ bins other than the "correct" one
static const double threshold = 0.08;

void
checkCQFreqOutput(const CQSpectrogram::RealBlock &output, double freq)
{
    
}

void
testCQFrequency(double freq)
{
    CQParameters params(sampleRate, cqmin, cqmax, bpo);
    CQSpectrogram cq(params, CQSpectrogram::InterpolateLinear);

    vector<double> input;
    for (int i = 0; i < duration; ++i) {
        input.push_back(sin((i * 2 * M_PI * freq) / sampleRate));
    }
    Window<double>(HanningWindow, duration).cut(input.data());

    CQSpectrogram::RealBlock output = cq.process(input);
    CQSpectrogram::RealBlock rest = cq.getRemainingOutput();
    output.insert(output.end(), rest.begin(), rest.end());

    checkCQFreqOutput(output, freq);
}

BOOST_AUTO_TEST_CASE(freq_0) { testCQFrequency(0); }
BOOST_AUTO_TEST_CASE(freq_5) { testCQFrequency(5); }
BOOST_AUTO_TEST_CASE(freq_10) { testCQFrequency(10); }
BOOST_AUTO_TEST_CASE(freq_15) { testCQFrequency(15); }
BOOST_AUTO_TEST_CASE(freq_20) { testCQFrequency(20); }
BOOST_AUTO_TEST_CASE(freq_25) { testCQFrequency(25); }
BOOST_AUTO_TEST_CASE(freq_30) { testCQFrequency(30); }
BOOST_AUTO_TEST_CASE(freq_35) { testCQFrequency(35); }
BOOST_AUTO_TEST_CASE(freq_40) { testCQFrequency(40); }
BOOST_AUTO_TEST_CASE(freq_45) { testCQFrequency(45); }
BOOST_AUTO_TEST_CASE(freq_50) { testCQFrequency(50); }

BOOST_AUTO_TEST_SUITE_END()

