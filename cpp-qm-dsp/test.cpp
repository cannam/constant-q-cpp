
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

