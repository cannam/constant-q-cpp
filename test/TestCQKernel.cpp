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

BOOST_AUTO_TEST_SUITE(TestCQKernel)

// Just some simple tests on kernel construction -- make sure it's the
// right size, etc

BOOST_AUTO_TEST_CASE(rate) {
    CQParameters params(123, 12, 65, 4);
    CQKernel k(params);
    BOOST_CHECK_EQUAL(k.getProperties().sampleRate, 123);
}

BOOST_AUTO_TEST_SUITE_END()

