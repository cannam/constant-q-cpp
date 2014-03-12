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

#include "ConstantQ.h"

#include <iostream>
#include <vector>

using std::vector;
using std::cerr;
using std::cout;
using std::endl;

#include <cstdio>

int main(int argc, char **argv)
{
    vector<double> in;

    for (int i = 0; i < 64; ++i) {
//	if (i == 0) in.push_back(1);
//	else in.push_back(0);
	in.push_back(sin(i * M_PI / 2.0));
    }

    ConstantQ k(8, 1, 4, 4);

    vector<vector<double> > out = k.process(in);
    vector<vector<double> > rest = k.getRemainingBlocks();

    out.insert(out.end(), rest.begin(), rest.end());

    cerr << "got " << out.size() << " back (" << out[0].size() << " in each?)" << endl;

    for (int b = 0; b < (int)out.size() / 8; ++b) {
	printf("\nColumns %d to %d:\n\n", b * 8, b * 8 + 7);
	for (int j = int(out[0].size()) - 1; j >= 0; --j) {
	    for (int i = 0; i < 8; ++i) {
		if (i + b * 8 < (int)out.size()) {
		    double v = out[i + b * 8][j];
		    if (v < 0.0001) printf("  0      ");
		    else printf("  %.4f ", out[i + b * 8][j]);
		}
	    }
	    printf("\n");
	}
    }
}

