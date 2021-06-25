#include <iostream>

#include "hf_wfn.hpp"
#include "hf_reader.hpp"

int main(int argc, char ** argv)
{
    hf_wfn master_wfn;
    read_2D_ints(argv[1],master_wfn,master_wfn.sints);
    read_2D_ints(argv[2],master_wfn,master_wfn.ke_ints);
    read_2D_ints(argv[3],master_wfn,master_wfn.eN_ints);
    read_4D_ints(argv[4],master_wfn);
    read_enuc(argv[5],master_wfn);

    master_wfn.make_core_H();

    master_wfn.orthogonalize_basis();
    master_wfn.make_initial_fock();
    master_wfn.diagonalize_fock();
    master_wfn.make_density_mat();
    master_wfn.evaluate_energy();
    
    //while(std::abs(master_wfn.total_e-master_wfn.last_e) > 0.000000001)
    while(!master_wfn.check_converged())
    {
        master_wfn.update_fock();
        master_wfn.diagonalize_fock();
        master_wfn.make_density_mat();
        master_wfn.evaluate_energy();        
    }
    
    master_wfn.print_mos();
    return 0;

}