#include <iostream>
#include <mpi.h>

int main(int argc, char** argv)
{

    int node;


    MPI_Init(&argc,&argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &node);

   std::cout << "Hello Wolrd from Node: " << node << std::endl;  
            
   MPI_Finalize();


    return 0;
}