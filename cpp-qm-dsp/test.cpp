
#include "ConstantQ.h"

#include <iostream>
#include <vector>

using std::vector;
using std::cerr;
using std::endl;

int main(int argc, char **argv)
{
    ConstantQ k(48000, 50, 24000, 24);

    vector<double> in(65536*4, 0.0);
    vector<vector<double> > out = k.process(in);


    cerr << "got " << out.size() << " back" << endl;
}

