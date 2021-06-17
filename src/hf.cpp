#include <iostream>

#include "hf_wfn.hpp"
#include "hf_reader.hpp"

int main(int argc, char ** argv)
{
    hf_wfn master_wfn;
    master_wfn.sints=read_2D_ints(argv[1],master_wfn);
    master_wfn.ke_ints=read_2D_ints(argv[2],master_wfn);
    master_wfn.eN_ints=read_2D_ints(argv[3],master_wfn);
    
    master_wfn.make_core_H();
    master_wfn.orthogonalize_basis();
    return 0;

}