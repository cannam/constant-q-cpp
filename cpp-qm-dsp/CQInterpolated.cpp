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

#include "CQInterpolated.h"

#include <iostream>
#include <stdexcept>

using std::vector;

using std::cerr;
using std::endl;

CQInterpolated::CQInterpolated(double sampleRate,
			       double minFreq, double maxFreq,
			       int binsPerOctave,
			       Interpolation interpolation) :
    m_cq(sampleRate, minFreq, maxFreq, binsPerOctave),
    m_interpolation(interpolation)
{
}

CQInterpolated::~CQInterpolated()
{
}

vector<vector<double> >
CQInterpolated::process(const vector<double> &td)
{
    return postProcess(m_cq.process(td), false);
}

vector<vector<double> >
CQInterpolated::getRemainingBlocks()
{
    return postProcess(m_cq.getRemainingBlocks(), true);
}

vector<vector<double> >
CQInterpolated::postProcess(const vector<vector<double> > &cq, bool insist)
{
    if (m_interpolation == None) {
	return cq;
    }

    int width = cq.size();

    for (int i = 0; i < width; ++i) {
	m_buffer.push_back(cq[i]);
    }
    
    if (m_interpolation == Hold) {
	return fetchHold(insist);
    } else {
	return fetchLinear(insist);
    }
}

vector<vector<double> >
CQInterpolated::fetchHold(bool)
{
    Grid out;
    
    int width = m_buffer.size();
    int height = getTotalBins();

    for (int i = 0; i < width; ++i) {
	
	vector<double> col = m_buffer[i];

	int thisHeight = col.size();
	int prevHeight = m_prevColumn.size();

	for (int j = thisHeight; j < height; ++j) {
	    if (j < prevHeight) {
		col.push_back(m_prevColumn[j]);
	    } else {
		col.push_back(0.0);
	    }
	}

	m_prevColumn = col;
	out.push_back(col);
    }

    m_buffer.clear();

    return out;
}

vector<vector<double> >
CQInterpolated::fetchLinear(bool insist)
{
    Grid out;

    //!!! This is surprisingly messy. I must be missing something.

    // We can only return any data when we have at least one column
    // that has the full height in the buffer, that is not the first
    // column.
    //
    // If the first col has full height, and there is another one
    // later that also does, then we can interpolate between those, up
    // to but not including the second full height column.  Then we
    // drop and return the columns we interpolated, leaving the second
    // full-height col as the first col in the buffer. And repeat as
    // long as enough columns are available.
    //
    // If the first col does not have full height, then (so long as
    // we're following the logic above) we must simply have not yet
    // reached the first full-height column in the CQ output, and we
    // can interpolate nothing.
    
    int width = m_buffer.size();
    int height = getTotalBins();

    if (width == 0) return out;
    
    int firstFullHeight = -1;
    int secondFullHeight = -1;

    for (int i = 0; i < width; ++i) {
	if ((int)m_buffer[i].size() == height) {
	    if (firstFullHeight == -1) {
		firstFullHeight = i;
	    } else if (secondFullHeight == -1) {
		secondFullHeight = i;
		break;
	    }
	}
    }

    if (firstFullHeight < 0) {
	if (insist) {
	    out = m_buffer;
	    m_buffer.clear();
	    return out;
	} else {
	    return out;
	}
    } else if (firstFullHeight > 0) {
	// can interpolate nothing, stash up to first full height & recurse
	out = Grid(m_buffer.begin(), m_buffer.begin() + firstFullHeight);
	m_buffer = Grid(m_buffer.begin() + firstFullHeight, m_buffer.end());
	Grid more = fetchLinear(insist);
	out.insert(out.end(), more.begin(), more.end());
	return out;
    } else if (secondFullHeight < 0) {
	// firstFullHeight == 0, but there is no second full height --
	// wait for it unless insist flag is set
	if (insist) {
	    out = m_buffer;
	    m_buffer.clear();
	    return out;
	} else {
	    return out;
	}
    } else {
	// firstFullHeight == 0 and secondFullHeight also valid. Can interpolate
	out = linearInterpolated(m_buffer, 0, secondFullHeight);
	m_buffer = Grid(m_buffer.begin() + secondFullHeight, m_buffer.end());
	Grid more = fetchLinear(insist);
	out.insert(out.end(), more.begin(), more.end());
	return out;
    }
}

vector<vector<double> >
CQInterpolated::linearInterpolated(const Grid &g, int x0, int x1)
{
    // g must be a grid with full-height columns at x0 and x1

    if (x0 >= x1) {
	throw std::logic_error("x0 >= x1");
    }
    if (x1 >= (int)g.size()) {
	throw std::logic_error("x1 >= g.size()");
    }
    if (g[x0].size() != g[x1].size()) {
	throw std::logic_error("x0 and x1 are not the same height");
    }

    int height = g[x0].size();
    int width = x1 - x0;

    Grid out(g.begin() + x0, g.begin() + x1);

    for (int y = 0; y < height; ++y) {

	int spacing = width;
	for (int i = 1; i < width; ++i) {
	    int thisHeight = g[x0 + i].size();
	    if (thisHeight > height) {
		throw std::logic_error("First column not full-height");
	    }
	    if (thisHeight > y) {
		spacing = i;
		break;
	    }
	}

	if (spacing < 2) continue;

	for (int i = 0; i + spacing <= width; i += spacing) {
	    for (int j = 1; j < spacing; ++j) {
		double proportion = double(j)/double(spacing);
		double interpolated = 
		    g[x0 + i][y] * (1.0 - proportion) +
		    g[x0 + i + spacing][y] * proportion;
		out[i + j].push_back(interpolated);
	    }
	}
    }

    return out;
}


	
