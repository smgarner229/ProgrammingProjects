#ifndef INPUTPARSER_HPP
#define INPUTPARSER_HPP

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "molecule.hpp"

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
        std::cout << "Error opening file: " << infile_name <<"\nExiting now\n";
        exit(1);
    }
    return;
}

#endif