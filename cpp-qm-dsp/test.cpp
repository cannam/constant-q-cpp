
#include "CQKernel.h"

#include <iostream>

int main(int argc, char **argv)
{
    CQKernel k(48000, 24000, 24);

    std::cerr << "Q = " << k.getProperties().Q << std::endl;

}

