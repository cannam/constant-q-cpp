/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/*
    Constant-Q library
    Copyright (c) 2013-2015 Queen Mary, University of London

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

#include "Chromagram.h"
#include "CQSpectrogram.h"
#include "Pitch.h"

#include <cstdio>

using namespace std;

Chromagram::Chromagram(Parameters params) :
    m_params(params),
    m_cq(0)
{
    int highestOctave = m_params.lowestOctave + m_params.octaveCount - 1;

    int midiPitchLimit = (1 + highestOctave) * 12 + 12; // C just beyond top
    double midiPitchLimitFreq = Pitch::getFrequencyForPitch
        (midiPitchLimit, 0, m_params.tuningFrequency);

    // Max frequency is frequency of the MIDI pitch just beyond the
    // top octave range (midiPitchLimit) minus one bin, then minus
    // floor(bins per semitone / 2)
    int bps = m_params.binsPerOctave / 12;
    m_maxFrequency = midiPitchLimitFreq /
        pow(2.0, (1.0 + floor(bps/2)) / m_params.binsPerOctave);

    // Min frequency is frequency of midiPitchLimit lowered by the
    // appropriate number of octaveCount.
    m_minFrequency = midiPitchLimitFreq /
        pow(2.0, m_params.octaveCount + 1);

    CQParameters p
        (params.sampleRate, m_minFrequency, m_maxFrequency, params.binsPerOctave);

    p.q = params.q;
    p.atomHopFactor = params.atomHopFactor;
    p.threshold = params.threshold;
    p.window = params.window;
    
    m_cq = new CQSpectrogram(p, CQSpectrogram::InterpolateLinear);
}

Chromagram::~Chromagram()
{
    delete m_cq;
}

bool
Chromagram::isValid() const
{
    return m_cq->isValid();
}

int
Chromagram::getColumnHop() const
{
    return m_cq->getColumnHop();
}

int
Chromagram::getLatency() const
{
    return m_cq->getLatency();
}

string
Chromagram::getBinName(int bin) const
{
    static const char *names[] = {
        "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"
    };

    float freq = m_cq->getBinFrequency(m_params.binsPerOctave - bin - 1);
    int note = Pitch::getPitchForFrequency(freq, 0, m_params.tuningFrequency);
    float nearestFreq =
        Pitch::getFrequencyForPitch(note, 0, m_params.tuningFrequency);
    
    char name[40];
    sprintf(name, "%d", bin);
    if (fabs(freq - nearestFreq) < 0.01) {
        return (name + std::string(" ") + names[note % 12]);
    } else {
        return (name);
    }
}

CQBase::RealBlock
Chromagram::process(const CQBase::RealSequence &data)
{
    return convert(m_cq->process(data));
}

CQBase::RealBlock
Chromagram::getRemainingOutput()
{
    return convert(m_cq->getRemainingOutput());
}

CQBase::RealBlock
Chromagram::convert(const CQBase::RealBlock &cqout)
{    
    CQBase::RealBlock chroma;

    int width = cqout.size();

    for (int i = 0; i < width; ++i) {

        CQBase::RealSequence column(m_params.binsPerOctave, 0.);

        // fold and invert to put low frequencies at the start

        int thisHeight = cqout[i].size();

	for (int j = 0; j < thisHeight; ++j) {
	    column[m_params.binsPerOctave - (j % m_params.binsPerOctave) - 1]
                += cqout[i][j];
	}

        chroma.push_back(column);
    }

    return chroma;
}

