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
    
    // Default constructor
    unit_vector(){};

    // Constructor for using two particles
    unit_vector(const particle A, const particle B, double magnitude)
    {
        x_dir=(A.x-B.x)/magnitude;
        y_dir=(A.y-B.y)/magnitude;
        z_dir=(A.z-B.z)/magnitude;
    };

    // Constructor for using three directions
    // Note this will Normalize the vector!
    unit_vector(const double x_dir, const double y_dir, const double z_dir):
    x_dir(x_dir), y_dir(y_dir), z_dir(z_dir) {};

    // Default destructor
    ~unit_vector(){};

};

struct vector
{
    struct unit_vector direction;
    double magnitude;

    // Constructor between two points
    vector(const particle & A, const particle & B)
    {
        magnitude = calc_distances(A,B);
        direction = unit_vector(A,B,magnitude);
    };

    // Constructor utilizing known magnitude and direction
    vector(const double magnitude, const unit_vector & direction):
    direction(direction),magnitude(magnitude) {};
};

double unit_dot_product(const unit_vector & ei, const unit_vector & ej)
{
    return ei.x_dir * ej.x_dir + ei.y_dir * ej.y_dir + ei.z_dir * ej.z_dir;
}

double vector_dot_product(const vector & ei, const vector & ej)
{
    return ei.magnitude * ej.magnitude * unit_dot_product(ei.direction,ej.direction);
}

unit_vector unit_cross_product(const unit_vector & ei, const unit_vector & ej)
{
    return unit_vector(ei.y_dir * ej.z_dir - ei.z_dir * ej.y_dir,
                       ei.z_dir * ej.x_dir - ei.x_dir * ej.z_dir,
                       ei.x_dir * ej.y_dir - ei.y_dir * ej.x_dir );
}

vector vector_cross_product(const vector & ei, const vector & ej)
{
    return vector(ei.magnitude * ej.magnitude, unit_cross_product(ei.direction,ej.direction));
}

double calc_angle(const particle & A, const particle & B, const particle & C)
{
    return std::acos(unit_dot_product(vector(A,B).direction,vector(C,B).direction)) * 180./M_PI;
}

//double calc_out_of_plane_angle(const particle & i, const particle & j, const particle & k, const particle & l)
//{
//    return std::asin(dot_product(cross_prodcut(unit_vector(j,k),unit_vector(l,k)),unit_vector(i,k))/(std::sin(std::acos(calc_angle(j,k,l)))))*180./M_PI;
//}

#endif