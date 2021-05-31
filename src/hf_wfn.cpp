#include "hf_wfn.hpp"
#include <cmath>

double * triangle_to_full_mat(double * triangle, const int & tri_size)
{
    double * full_mat;
    int full_size = std::pow((-1+std::pow(1+8*tri_size,0.5))/2,1);
    full_mat = new double[full_size];
    size_t counteri=0,counterj=0;

    for(size_t i = 0; i < tri_size; i++)
    {
        if (counteri == counterj)
        {
            full_mat[counteri*full_size+counterj]=triangle[i];
            counteri++;
            counterj=0;
        }
        else
        {
            full_mat[counteri*full_size+counterj]=triangle[i];
            full_mat[counterj*full_size+counteri]=triangle[i];
            counterj++;
        }
    }

    for(size_t i = 0; i < full_size; i++)
    {
        std::cout << '\n';
        for(size_t j = 0; j < full_size; j++)
        {
            std::cout << '\t' << full_mat[i*full_size+j];
        }
    }
    //delete [] triangle;
    return full_mat;
}

double two_center_integral::operator()(int i, int j)
{

}

bool two_center_integral::operator==(const two_center_integral & other)
{

}

void hf_wfn::print_sints()
{
    for (size_t i = 0; i < mat_size; i++)
    {
        std::cout << sints[i] << std::endl;
    }
}

void hf_wfn::make_core_H()
{
    core_H = new double[mat_size];
    for(size_t j = 0; j < mat_size; j++)
    {
        core_H[j] = ke_ints[j] + eN_ints[j];
    }
    for(size_t i = 0; i < mat_size; i++)
    {
        std::cout << core_H[i] << std::endl;
    }
    return;
}

void hf_wfn::orthogonalize_basis()
{
    char charv='V';
    char charL='L';
    int N = (-1+std::pow(1+8*mat_size,0.5))/2;
    
}