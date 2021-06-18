#include "hf_wfn.hpp"
#include <cmath>
#include <iomanip>

#define FIELDWIDTH 10

extern "C" {
    extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
    extern int dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
    extern int dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
}

static void print_mat(double * & mat, const int & stride)
{
    for(size_t i = 0; i < stride; i++)
    {
        std::cout << '\n';
        for(size_t j = 0; j < stride; j++)
        {
            if(std::abs(mat[i*stride+j])>1e-9) std::cout << std::setw(FIELDWIDTH) << std::setprecision(5) << mat[i*stride+j];
            else std::cout << std::setw(FIELDWIDTH) << std::setprecision(5) << 0.0;
        }
    }
    std::cout << "\n\n";
    return;
}

void print_triangle_as_full_mat(double * & triangle, const int & tri_size)
{
    double * temp = nullptr;
    size_t full_size = std::pow((-1+std::pow(1+8*tri_size,0.5))/2,1);
    temp=triangle_to_full_mat(triangle,tri_size);

    print_mat(temp,full_size);

    delete [] temp;
    return;
}

double * triangle_to_full_mat(double * & triangle, const int & tri_size)
{
    double * full_mat = nullptr;
    int full_size = std::pow((-1+std::pow(1+8*tri_size,0.5))/2,1);
    full_mat = new double[full_size*full_size];
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

    if(false)
    {
        for(size_t i = 0; i < full_size; i++)
        {
            std::cout << '\n';
            for(size_t j = 0; j < full_size; j++)
            {
                std::cout << '\t' << full_mat[i*full_size+j];
            }
        }
    }
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
    if(false)
    {
    for(size_t i = 0; i < mat_size; i++)
    {
        std::cout << core_H[i] << std::endl;
    }
    }
    return;
}

void hf_wfn::orthogonalize_basis()
{

    int full_size = std::pow((-1+std::pow(1+8*mat_size,0.5))/2,1);
    double * temp = triangle_to_full_mat(sints,mat_size);
    char N = 'N';
    char V = 'V';
    int one = 1;
    double * eigRE = new double[full_size]{0.0};
    double * eigIM = new double[full_size]{0.0};
    double * eigL = new double[full_size*full_size]{0.0};
    double * eigR = new double[full_size*full_size]{0.0};
    double * work = new double[full_size*full_size]{0.0};
    int lwork = 4 * full_size;
    int info = 0;
    dgeev_(&N,&V,&full_size,temp,&full_size,eigRE,eigIM,eigL,&full_size,eigR,&full_size,work,&lwork,&info);
    std::cout << "info: " << info << std::endl;
    print_mat(eigR,full_size);

    double * eigval_mat = new double[full_size*full_size]{0.0};
    for (size_t i = 0; i < full_size; i++)
    {
        eigval_mat[i*full_size+i] = std::pow(eigRE[i],-0.5);
        //eigval_mat[i*full_size+i] = 1.0;
        std::cout << eigRE[i] << std::endl;
    }
    print_mat(eigval_mat,full_size);
    print_mat(eigR,full_size);
    
    double * orthmat = new double[full_size*full_size]{0.0};
    char T = 'T';
    double oned =1.;
    double zerod =0.;
    double * final_mat  = new double[full_size*full_size]{0.0};

    dgemm_(&N,&N,&full_size,&full_size,&full_size,&oned,eigR,&full_size,eigval_mat,&full_size,&oned,orthmat,&full_size);
    print_mat(orthmat,full_size);
    dgemm_(&N,&T,&full_size,&full_size,&full_size,&oned,orthmat,&full_size,eigR,&full_size,&oned,final_mat,&full_size);
    print_mat(final_mat,full_size);
    

    delete[] temp;
    
}