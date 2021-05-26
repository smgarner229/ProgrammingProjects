#ifndef INPUTPARSER_HPP
#define INPUTPARSER_HPP

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include "molecule.hpp"

static void file_open_error(char * file_name)
{
    std::cout << "Cannot open file: " << file_name << "\nExitting now!\n";
    exit(1);
}

void parse_geom(char* infile_name, molecule & mol)
{
    double read_val;

    std::vector<double> read_pos;
    std::ifstream infile(infile_name);

    int natoms;

    if(infile)
    {
        infile >> natoms;

        while (true)
        {
            infile >> read_val;
            if(infile.eof())
            { 
                mol.add_neucleus(read_pos[0],read_pos[1],read_pos[2],read_pos[3]);
                break; 
            }
            if (read_pos.size()==4)
            {
                mol.add_neucleus(read_pos[0],read_pos[1],read_pos[2],read_pos[3]);
                read_pos.clear();
                read_pos.push_back(read_val);
            }
            else
            {
                read_pos.push_back(read_val);
            }
        }
    }
    else
    {
        file_open_error(infile_name);
    }
    return;
}

void parse_hessian(char * infile_name, molecule & mol)
{
    int read_natoms;
    int counter = 0;
    std::ifstream infile(infile_name);

    if (infile) // File opened successfully
    {
        infile >> read_natoms;
        int total_size = (int)std::pow(3*read_natoms,2.0);
        mol.hessian = new double[total_size];
        while(counter < total_size)
        {
            infile >> mol.hessian[counter];
            std::cout << mol.hessian[counter] << std::endl;
            counter++;
        }
    }
    else
    {
        file_open_error(infile_name);
    }

    return;
}

#endif