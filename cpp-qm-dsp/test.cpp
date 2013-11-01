
#include "CQKernel.h"

#include <iostream>

int main(int argc, char **argv)
{
    CQKernel k(96000, 48000, 24);

    std::cerr << "Q = " << k.getProperties().Q << std::endl;

}

