/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "cq/CQSpectrogram.h"

#include "dsp/Window.h"

#include <cmath>
#include <vector>
#include <iostream>

using std::vector;
using std::cerr;
using std::endl;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(TestCQTime)

// Principle: Run a Dirac impulse through the CQ transform and check
// that its output has all the peak bins aligned correctly in time.

// Set up fs/2 = 50, frequency range 10 -> 40 i.e. 2 octaves, fixed
// duration of 2 seconds
static const double sampleRate = 100;
static const double cqmin = 10;
static const double cqmax = 40;
static const double bpo = 4;
static const int duration = sampleRate * 2;

// Threshold below which to ignore a column completely
static const double threshold = 0.08;

void
testCQTime(double t)
{
    vector<CQSpectrogram::Interpolation> interpolationTypes;
    interpolationTypes.push_back(CQSpectrogram::InterpolateZeros);
    interpolationTypes.push_back(CQSpectrogram::InterpolateHold);
    interpolationTypes.push_back(CQSpectrogram::InterpolateLinear);

    for (int k = 0; k < int(interpolationTypes.size()); ++k) {

        CQParameters params(sampleRate, cqmin, cqmax, bpo);
        CQSpectrogram cq(params, interpolationTypes[k]);

        BOOST_CHECK_EQUAL(cq.getBinsPerOctave(), bpo);
        BOOST_CHECK_EQUAL(cq.getOctaves(), 2);

        //!!! generate input signal
        vector<double> input;


        CQSpectrogram::RealBlock output = cq.process(input);
        CQSpectrogram::RealBlock rest = cq.getRemainingOutput();
        output.insert(output.end(), rest.begin(), rest.end());

        BOOST_CHECK_EQUAL(output[0].size(), 
                          cq.getBinsPerOctave() * cq.getOctaves());

        //!!! test output signal
    }
}

BOOST_AUTO_TEST_CASE(time_zero) { testCQTime(0); }
BOOST_AUTO_TEST_CASE(time_half) { testCQTime(0.5); }
BOOST_AUTO_TEST_CASE(time_one) { testCQTime(1.0); }
BOOST_AUTO_TEST_CASE(time_two) { testCQTime(2.0); }

BOOST_AUTO_TEST_SUITE_END()

