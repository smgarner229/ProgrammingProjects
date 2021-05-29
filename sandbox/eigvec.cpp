#include <iostream>
#include <cmath>

extern "C" {
    extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
    extern int dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
    }

void print_mat(double * printme);

int main(int argc, char ** argv)
{
    double * matrix;
    matrix = new double[9] {1.0,2.0,0.0,
                            2.0,2.0,0.0,
                            0.0,0.0,3.0};

    char N = 'N';
    char V = 'V';
    int n = 3;
    double * eigReal = new double[n];
    double * eigImag = new double[n];
    double * lefteigVec = nullptr;
    double * righteigVec = new double[n*n]{0.0};
    int lwork = 6*n;
    double * work = new double[lwork];
    int info = 0;

    print_mat(matrix);

    dgeev_(&N,&V,&n,matrix,&n,eigReal,eigImag,lefteigVec,&n,righteigVec,&n,work,&lwork,&info);

    print_mat(matrix);
    print_mat(righteigVec);

    for (size_t j = 0; j < 3; j++)
    {
        std::cout << eigReal[j] << "\n";
    }
    std::cout << "\n";

    double * righttrans = new double[n*n]{0.0};
    double * store_solution = new double[n*n]{0.0};
    for(size_t i = 0; i < 9; i++)
    {
        righttrans[i] = righteigVec[i];
    }
    
    char TA = 'N';
    char TB = 'T';
    double one = 1.;
    double zero = 0.;

    print_mat(righteigVec);
    print_mat(righttrans);
    
    dgemm_(&TA,&TB,&n,&n,&n,&one,righteigVec,&n,righttrans,&n,&zero,store_solution,&n);
    print_mat(store_solution);


    delete work;
    delete eigReal;
    delete eigImag;
    delete righttrans;
    delete righteigVec;
    delete matrix;
    delete store_solution;

    return 0;
}

void print_mat(double * printme)
{
    for (size_t i = 0; i < 9; i++)
    {
        if (i%3==0) std::cout << "\n";
        std::cout << printme[i] << "\t";
    }
    std::cout << "\n";
    return;
}