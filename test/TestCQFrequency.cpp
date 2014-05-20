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

BOOST_AUTO_TEST_SUITE(TestCQFrequency)

// The principle here is to feed a single windowed sinusoid into a
// small CQ transform and check that the output has its peak bin at
// the correct frequency. 

// Set up fs/2 = 50, frequency range 10 -> 40 i.e. 2 octaves, fixed
// duration of 2 seconds
static const double sampleRate = 100;
static const double cqmin = 11.8921;
static const double cqmax = 40;
static const double bpo = 4;
static const int duration = sampleRate * 2;

// Threshold below which to ignore a column completely
static const double threshold = 0.08;

int
binForFrequency(double freq)
{
    int bin = (bpo * 2) - round(bpo * log2(freq / cqmin)) - 1;
    return bin;
}

void
checkCQFreqColumn(int i, vector<double> column,
                  double freq, CQSpectrogram::Interpolation interp)
{
    double maxval = 0.0;
    int maxidx = -1;
    int height = column.size();

    int nonZeroHeight = ((i % 2 == 1) ? height/2 : height);

    for (int j = 0; j < nonZeroHeight; ++j) {
        if (j == 0 || column[j] > maxval) {
            maxval = column[j];
            maxidx = j;
        }
    }

    int expected = binForFrequency(freq);
    if (maxval < threshold) {
        return; // ignore these columns at start and end
    } else if (expected < nonZeroHeight && maxidx != expected) {
        cerr << "ERROR: In column " << i << " with interpolation " << interp
             << ", maximum value for frequency " << freq
             << "\n       found at index " << maxidx
             << " (expected index " << expected << ")" << endl;
        cerr << "column contains: ";
        for (int j = 0; j < height; ++j) {
            cerr << column[j] << " ";
        }
        cerr << endl;
        BOOST_CHECK_EQUAL(maxidx, expected);
    }
}

void
testCQFrequency(double freq)
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
        BOOST_CHECK_CLOSE(cq.getBinFrequency(0), 40, 1e-10);
        BOOST_CHECK_CLOSE(cq.getBinFrequency(4), 20, 1e-10);
        BOOST_CHECK_CLOSE(cq.getBinFrequency(7), cqmin, 1e-3);

        vector<double> input;
        for (int i = 0; i < duration; ++i) {
            input.push_back(sin((i * 2 * M_PI * freq) / sampleRate));
        }
        Window<double>(HanningWindow, duration).cut(input.data());

        CQSpectrogram::RealBlock output = cq.process(input);
        CQSpectrogram::RealBlock rest = cq.getRemainingOutput();
        output.insert(output.end(), rest.begin(), rest.end());

        BOOST_CHECK_EQUAL(output[0].size(), 
                          cq.getBinsPerOctave() * cq.getOctaves());

        for (int i = 0; i < int(output.size()); ++i) {
            checkCQFreqColumn(i, output[i], freq, interp);
        }
    }
}

BOOST_AUTO_TEST_CASE(freq_11) { testCQFrequency(11); }
BOOST_AUTO_TEST_CASE(freq_17) { testCQFrequency(17); }
BOOST_AUTO_TEST_CASE(freq_24) { testCQFrequency(24); }
BOOST_AUTO_TEST_CASE(freq_27) { testCQFrequency(27); }
BOOST_AUTO_TEST_CASE(freq_33) { testCQFrequency(33); }
BOOST_AUTO_TEST_CASE(freq_40) { testCQFrequency(40); }

BOOST_AUTO_TEST_SUITE_END()

