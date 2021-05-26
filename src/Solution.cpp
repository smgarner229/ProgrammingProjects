#include <iostream>
#include <fstream>
#include "inputparser.hpp"
#include "molecule.hpp"

#define MAXLINE 80

int main(int argc, char** argv)
{
    molecule mastermol;
    parse_geom(argv[1],mastermol);
    parse_hessian(argv[2],mastermol);
    std::cout << mastermol;
    mastermol.calc_bond_lengths();
    mastermol.calc_bond_angles();
    mastermol.calc_outofplane_angle();
    mastermol.calc_torsion_angle();
    mastermol.calc_center_of_mass();
    mastermol.calc_inertial_tensor();
    mastermol.vibrational_analysis();

    return 0;
}