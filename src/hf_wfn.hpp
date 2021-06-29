#ifndef HF_WFN_HPP
#define HF_WFN_HPP

#include <map>
#include <utility>
#include <iostream>

#include "tei_handler.hpp"

void print_triangle_as_full_mat(double * & triangle, const int & tri_size);
double * triangle_to_full_mat(double * & triangle, const int & tri_size);

class hf_wfn
{
    public:
        size_t mat_size = -1;
        int nao;
        int niter=0;
        double enuc;
        double * sints = nullptr;
        double * ke_ints = nullptr;
        double * eN_ints = nullptr;
        double * core_H = nullptr;
        double * sym_orth_mat = nullptr;
        two_electron_integral_handler teis;
        double * fock = nullptr;
        double * c_mat = nullptr;
        double * density_mat = nullptr;
        double * last_density_mat = nullptr;
        double total_e;
        double * inbasis_fock = nullptr;
        double last_e = 0.0;
        double * orb_energies;

    hf_wfn(){};
    ~hf_wfn(){delete[] sints; 
              delete[] ke_ints; 
              delete[] eN_ints; 
              delete[] core_H; 
              delete[] sym_orth_mat; 
              delete[] fock;
              delete[] c_mat; 
              delete[] density_mat;};

    void print_sints();
    void make_core_H();
    void orthogonalize_basis();
    void make_initial_fock();
    void make_density_mat();
    void build_fock();
    void diagonalize_fock();
    void evaluate_energy();
    void update_fock();
    void print_mos();
    void fock_in_mo_basis();
    bool check_converged();

};

#endif