#ifndef HF_WFN_HPP
#define HF_WFN_HPP

#include <map>
#include <utility>
#include <iostream>

void print_triangle_as_full_mat(double * & triangle, const int & tri_size);
double * triangle_to_full_mat(double * & triangle, const int & tri_size);

class two_center_integral
{
    private:
        int i,j;
        double integral;
    public:
        double operator ()(int i, int j);
        bool operator ==(const two_center_integral & other);
};

class overlap_integral : public two_center_integral
{
    public:
        double operator()(int i, int j);
};

class kinetic_integral : public two_center_integral
{
    public:
        double operator()(int i, int j);
};



class hf_wfn
{
    public:
        size_t mat_size = -1;
        double enuc;
        double * sints = nullptr;
        double * ke_ints = nullptr;
        double * eN_ints = nullptr;
        double * core_H = nullptr;
    hf_wfn(){};
    ~hf_wfn(){delete [] sints; delete [] ke_ints; delete [] eN_ints; delete [] core_H;};

    void print_sints();
    void make_core_H();
    void orthogonalize_basis();
};

#endif