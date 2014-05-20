/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "cq/CQKernel.h"

#include <cmath>
#include <vector>
#include <iostream>

using std::vector;
using std::cerr;
using std::endl;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

static int rate = 123;
static int max = 60;
static int min = 12;
static int bpo = 4;

BOOST_AUTO_TEST_SUITE(TestCQKernel)

// Just some simple tests on kernel construction -- make sure it's the
// right size, etc

BOOST_AUTO_TEST_CASE(sampleRate) {
    CQParameters params(rate, min, max, bpo);
    CQKernel k(params);
    BOOST_CHECK_EQUAL(k.getProperties().sampleRate, rate);
}

BOOST_AUTO_TEST_CASE(binsPerOctave) {
    CQParameters params(rate, min, max, bpo);
    CQKernel k(params);
    BOOST_CHECK_EQUAL(k.getProperties().binsPerOctave, bpo);
}

BOOST_AUTO_TEST_CASE(maxFrequency) {
    CQParameters params(rate, min, max, bpo);
    CQKernel k(params);
    BOOST_CHECK_EQUAL(k.getProperties().maxFrequency, max);
}

BOOST_AUTO_TEST_CASE(minFrequency) {
    CQParameters params(rate, min, max, bpo);
    CQKernel k(params);
    BOOST_CHECK_CLOSE(k.getProperties().minFrequency,
		      (max / 2.0) * pow(2, 1.0/bpo),
		      1e-8);
}

BOOST_AUTO_TEST_CASE(atomsPerFrame) {
    CQParameters params(rate, min, max, bpo);
    CQKernel k(params);
    BOOST_CHECK_EQUAL(k.getProperties().atomsPerFrame, 5);
}

BOOST_AUTO_TEST_CASE(fftSize) {
    CQParameters params(rate, min, max, bpo);
    CQKernel k(params);
    BOOST_CHECK_EQUAL(k.getProperties().fftSize, 32);
}

BOOST_AUTO_TEST_SUITE_END()

