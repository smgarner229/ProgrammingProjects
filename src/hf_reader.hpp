#include <iostream>
#include <fstream>

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

double * read_2D_ints(const char * infile_name, hf_wfn & wfn)
{
    double * store_type;
    if (wfn.mat_size == -1)
    {
        wfn.mat_size = get_mat_size(infile_name);
    }
    store_type = new double[wfn.mat_size];
    size_t index = 0;
    std::ifstream open_file = file_opener(infile_name);
    int temp;
    while(index < wfn.mat_size)
    {
        open_file >> temp;
        open_file >> temp;
        open_file >> store_type[index++];
    }
    open_file.close();
    for (size_t j = 0; j < wfn.mat_size; j++)
    {
        std::cout << store_type[j] << std::endl;
    }
    return store_type;
}
