#include <iostream>
#include <fstream>
#include <cmath>

#include "hf_wfn.hpp"

std::ifstream file_opener(const char * infile_name)
{
    std::ifstream open_file(infile_name);
    if(! open_file)
    {
        std::cout << "Error opening file: " << infile_name << "\nExitting now\n";
        exit(1);
    }
    std::cout << "Opened file: " << infile_name << "\n";
    return open_file;
}

void read_enuc(const char * infile_name, hf_wfn & wfn)
{
    std::ifstream open_file = file_opener(infile_name);
    double enuc;
    open_file >> enuc;
    wfn.enuc = enuc;
    open_file.close();
    return;
}

int get_mat_size(const char * infile_name)
{
    int findex, sindex;
    double integral;
    std::ifstream open_file = file_opener(infile_name);
    while(!open_file.eof())
    {
        open_file >> findex;
        open_file >> sindex;
        open_file >> integral;
    }
    open_file.close();
    return findex * (findex + 1) / 2;
}

void read_2D_ints(const char * infile_name, hf_wfn & wfn, double * & storage_location)
{
    if (wfn.mat_size == -1)
    {
        wfn.mat_size = get_mat_size(infile_name);
    }
    storage_location = new double[wfn.mat_size];

    size_t index = 0;
    std::ifstream open_file = file_opener(infile_name);
    int temp;
    while(index < wfn.mat_size)
    {
        open_file >> temp;
        open_file >> temp;
        open_file >> storage_location[index++];
    }
    open_file.close();

    if(false)
    {
        for (size_t j = 0; j < wfn.mat_size; j++)
        {
            std::cout << storage_location[j] << std::endl;
        }
    }
    return;

    // At this point, we're stored into a 1D array for the lower triangle.
    // In the future just use this, but for now let's use the larger 2D
    // Full version since these matrices are small enough
    /*
    if(false){
    double * full_mat = NULL;
    int full_size = (-1+std::pow(1+8*wfn.mat_size,0.5))/2;

    full_mat = new double[full_size*full_size];
    size_t counteri=0,counterj=0;

    for(size_t i = 0; i < wfn.mat_size; i++)
    {
        std::cout << "i is\t" << i << std::endl;
        if (counteri == counterj)
        {
            full_mat[counteri*full_size+counterj]=storage_location[i];
            std::cout << "SpotD: " << counteri*full_size + counterj << "\t" << counterj*full_size + counteri << "\n";
            counteri++;
            counterj=0;
        }
        else
        {
            std::cout << "Spot: " << counteri*full_size + counterj << "\t" << counterj*full_size + counteri << "\n";
            full_mat[counteri*full_size+counterj]=storage_location[i];
            full_mat[counterj*full_size+counteri]=storage_location[i];
            counterj++;
        }
        std::cout << counteri << "\t" << counterj << "\n"; 
    }
    for(size_t i = 0; i < full_size; i++)
    {
        std::cout << '\n';
        for(size_t j = 0; j < full_size; j++)
        {
            std::cout << '\t' << full_mat[i*full_size+j];
        }

        }    
    }
    */
    //return full_mat;
}
