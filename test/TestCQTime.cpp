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

        CQSpectrogram::Interpolation interp = interpolationTypes[k];

        CQParameters params(sampleRate, cqmin, cqmax, bpo);
        CQSpectrogram cq(params, interp);

        BOOST_CHECK_EQUAL(cq.getBinsPerOctave(), bpo);
        BOOST_CHECK_EQUAL(cq.getOctaves(), 2);

        //!!! generate input signal
        vector<double> input(duration, 0.0);
        int ix = int(floor(t * sampleRate));
        if (ix >= duration) ix = duration-1;
        input[ix] = 1.0;

        CQSpectrogram::RealBlock output = cq.process(input);
        CQSpectrogram::RealBlock rest = cq.getRemainingOutput();
        output.insert(output.end(), rest.begin(), rest.end());

        BOOST_CHECK_EQUAL(output[0].size(), 
                          cq.getBinsPerOctave() * cq.getOctaves());

        vector<int> peaks;
        double eps = 1e-8;

        for (int j = 0; j < int(output[0].size()); ++j) {

            int maxidx = -1;
            double max = 0.0;
            for (int i = 0; i < int(output.size()); ++i) {
                double value = output[i][j];
                if (i == 0 || value + eps > max) {
                    max = value;
                    maxidx = i;
                }
            }

            peaks.push_back(maxidx);
        }

        for (int j = 1; j < int(peaks.size()); ++j) {
            int oct = j / bpo;
            int spacing = (1 << oct);
            int actual = peaks[j]/spacing;
            int expected = int(round(double(peaks[0])/spacing));
            if (actual != expected) {
                cerr << "ERROR: In row " << j << " (bin freq "
                     << cq.getBinFrequency(j) << "), interpolation " << interp
                     << ", maximum value for time " << t 
                     << "\n       found at index " << peaks[j] 
                     << " of " << output.size() << " which does not align with"
                     << " highest frequency\n       bin peak at " << peaks[0]
                     << " given octave spacing of " << spacing
                     << "\n       [latency = " << cq.getLatency()
                     << ", hop = " << cq.getColumnHop() << ", duration = "
                     << duration << ", ix = " << ix << "]" << endl;
                cerr << "row contains: ";
                for (int i = 0; i < int(output.size()); ++i) {
                    if (i == expected * spacing) cerr << "*";
                    if (i == peaks[j]) cerr << "**";
                    cerr << output[i][j] << " ";
                }
                cerr << endl;

                BOOST_CHECK_EQUAL(actual, expected);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(time_zero) { testCQTime(0); }
BOOST_AUTO_TEST_CASE(time_half) { testCQTime(0.5); }
BOOST_AUTO_TEST_CASE(time_one) { testCQTime(1.0); }
BOOST_AUTO_TEST_CASE(time_two) { testCQTime(2.0); }

BOOST_AUTO_TEST_SUITE_END()

