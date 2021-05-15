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
    mastermol.calc_bond_lengths();
    mastermol.calc_bond_angles();

    return 0;
}