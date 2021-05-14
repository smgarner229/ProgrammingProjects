#include <iostream>
#include <fstream>
#include "inputparser.hpp"
#include "molecule.hpp"

#define MAXLINE 80

int main(int argc, char** argv)
{
    molecule mastermol;
    parse_geom(argv[1],mastermol);
    std::cout << mastermol;

    return 0;
}