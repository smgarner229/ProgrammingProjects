#ifndef GEOMETRY_OPPS_HPP
#define GEOMETRY_OPP_HPP

#include <cmath>
#include <vector>

#include "molecule.hpp"

double calc_distances(const particle A, const particle B)
{
    return std::sqrt(std::pow(A.x-B.x,2.)+std::pow(A.y-B.y,2.)+std::pow(A.z-B.z,2.));
}

struct unit_vector
{
    double x_dir, y_dir, z_dir;
    
    // Constructor for using two particles
    unit_vector(const particle A, const particle B)
    {
        double magnitue = calc_distances(A,B);
        x_dir=-1.0*(B.x-A.x)/magnitue;
        y_dir=-1.0*(B.y-A.y)/magnitue;
        z_dir=-1.0*(B.z-A.z)/magnitue;
    };

    // Constructor which takes the dot prodcut.  
    // Probably bad practice to not have this explicity in the function call name
    unit_vector(const unit_vector &ei, const unit_vector & ej)
    {
        x_dir=ei.y_dir*ej.z_dir-ei.z_dir*ej.y_dir;
        y_dir=ei.z_dir*ej.x_dir-ei.x_dir*ej.z_dir;
        z_dir=ei.x_dir*ej.y_dir-ei.y_dir*ej.x_dir;
    };
};

double dot_product(const unit_vector & ei, const unit_vector & ej)
{
    return ei.x_dir * ej.x_dir + ei.y_dir * ej.y_dir + ei.z_dir * ej.z_dir;
}

unit_vector cross_prodcut(const unit_vector & ei, const unit_vector & ej)
{
    return unit_vector(ei,ej);
}


double calc_angle(const particle & B, const particle & A, const particle & C)
{
    return std::acos(dot_product(unit_vector(A,B), unit_vector(C,B)))*180./M_PI;
}

double calc_out_of_plane_angle(const particle & i, const particle & j, const particle & k, const particle & l)
{
    return std::asin(dot_product(cross_prodcut(unit_vector(j,k),unit_vector(l,k)),unit_vector(i,k))/(std::sin(std::acos(calc_angle(j,k,l)))))*180./M_PI;
}

#endif