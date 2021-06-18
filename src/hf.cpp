#include <iostream>

#include "hf_wfn.hpp"
#include "hf_reader.hpp"

int main(int argc, char ** argv)
{
    hf_wfn master_wfn;
    read_2D_ints(argv[1],master_wfn,master_wfn.sints);
    print_triangle_as_full_mat(master_wfn.sints,master_wfn.mat_size);
    read_2D_ints(argv[2],master_wfn,master_wfn.ke_ints);
    print_triangle_as_full_mat(master_wfn.ke_ints,master_wfn.mat_size);
    read_2D_ints(argv[3],master_wfn,master_wfn.eN_ints);
    print_triangle_as_full_mat(master_wfn.eN_ints,master_wfn.mat_size);

    master_wfn.make_core_H();

    read_4D_ints(argv[4],master_wfn);

    master_wfn.orthogonalize_basis();
    
    return 0;

}