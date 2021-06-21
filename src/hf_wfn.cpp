#include "hf_wfn.hpp"
#include <cmath>
#include <iomanip>

#define FIELDWIDTH 10

extern "C" {
    extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
    extern int sgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
    extern int dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
    extern int dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
}

static void swap(double * one, double * two)
{
    double temp = *one;
    *one = *two;
    *two = temp;
}

static void print_mat(double * & mat, const int & stride, bool trans = false)
{
    size_t index;
    for(size_t i = 0; i < stride; i++)
    {
        std::cout << '\n';
        for(size_t j = 0; j < stride; j++)
        {
            trans ? index = j * stride + i : index = i * stride + j;
            if(std::abs(mat[index])>1e-9) std::cout << std::setw(FIELDWIDTH) << std::setprecision(5) << mat[index];
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
    return 0.0;
}

bool two_center_integral::operator==(const two_center_integral & other)
{
    return false;
}

void hf_wfn::print_sints()
{
    for (size_t i = 0; i < mat_size; i++)
    {
        std::cout << sints[i] << std::endl;
    }
    return;
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
    int full_size = std::pow((-1+std::pow(1+8*mat_size,0.5))/2,1);
    core_H = triangle_to_full_mat(core_H,mat_size);
    return;
}

void hf_wfn::orthogonalize_basis()
{
    int full_size = std::pow((-1+std::pow(1+8*mat_size,0.5))/2,1);
    double * temp = new double[full_size*full_size];
    temp = triangle_to_full_mat(sints,mat_size);
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

    // Diagonalize the overlap matrix
    // Eigenvalues are stored in eigRE (since we know these should be real.....for now....)
    // Associated Eigenvectors are stored in eigR
    dgeev_(&N,&V,&full_size,temp,&full_size,eigRE,eigIM,eigL,&full_size,eigR,&full_size,work,&lwork,&info);

    double * eigval_mat = new double[full_size*full_size]{0.0};
    for (size_t i = 0; i < full_size; i++)
    {
        eigval_mat[i*full_size+i] = std::pow(eigRE[i],-0.5);
    }

    double * orthmat = new double[full_size*full_size]{0.0};
    char T = 'T';
    double oned =1.;
    double zerod =0.;

    // Reuse the temp matrix to store the result of L*S
    dgemm_(&N,&N,&full_size,&full_size,&full_size,&oned,eigR,&full_size,eigval_mat,&full_size,&zerod,temp,&full_size);
    sym_orth_mat = new double[full_size*full_size]{0.0};

    // Carry out (L*S)*L^T, store the result in sym_orth_mat
    dgemm_(&N,&T,&full_size,&full_size,&full_size,&oned,temp,&full_size,eigR,&full_size,&zerod,sym_orth_mat,&full_size);

    //delete[] orthmat;
    delete[] eigRE;
    delete[] eigIM;
    delete[] eigL;
    delete[] eigR;
    delete[] work;
    delete[] eigval_mat;

    delete[] temp;
    return;
}

void hf_wfn::make_initial_fock()
{
    int full_size = std::pow((-1+std::pow(1+8*mat_size,0.5))/2,1);
    fock= new double[full_size*full_size]{0.0};
    double * temp = new double[full_size*full_size]{0.0};

    char N = 'N';
    char T = 'T';
    double oned =1.;
    double zerod =0.;

    inbasis_fock = new double[full_size*full_size]{0.0};
    for(size_t i = 0; i < full_size*full_size; i++)
        inbasis_fock[i]=core_H[i];

    dgemm_(&T,&N,&full_size,&full_size,&full_size,&oned,sym_orth_mat,&full_size,core_H,&full_size,&zerod,temp,&full_size);
    dgemm_(&N,&N,&full_size,&full_size,&full_size,&oned,temp,&full_size,sym_orth_mat,&full_size,&zerod,fock,&full_size);

    delete[] temp;
    return;
}

static void sort_mos(int nmo, double * & mo_energy, double * & mos)
{
    double * sorted_mos = new double[nmo*nmo]{0.0};
    double * sorted_mo_e = new double[nmo]{0.0};
    double tempmo;

    int min_idx;
    size_t j;
    for(size_t i = 0; i < nmo-1; i++)
    {
        min_idx=i;
        for(j = i+1; j < nmo; j++)
        {
            if(mo_energy[j]<mo_energy[min_idx])
            {
                min_idx=j;
            }
        }
        if(min_idx!=i)
        {
            swap(&mo_energy[min_idx],&mo_energy[i]);
            for(size_t k = 0; k < nmo; k++)
            {
                tempmo=mos[min_idx*nmo + k];
                mos[min_idx*nmo+k]=mos[i*nmo+k];
                mos[i*nmo+k]=tempmo;
            }
        }
    }
    return;
}

void hf_wfn::diagonalize_fock()
{
    int full_size = std::pow((-1+std::pow(1+8*mat_size,0.5))/2,1);
    double * temp = new double[full_size*full_size];
    temp = triangle_to_full_mat(sints,mat_size);
    char N = 'N';
    char V = 'V';
    int one = 1;
    double * eigRE = new double[full_size]{0.0};
    double * eigIM = new double[full_size]{0.0};
    double * eigL = new double[full_size*full_size]{0.0};
    double * eigR = new double[full_size*full_size]{0.0};
    double * work = new double[full_size*full_size]{0.0};
    int lwork = 6 * full_size;
    int info = 0;
    c_mat = new double[full_size*full_size];
    double * copy_foc = new double[full_size*full_size];
    for(size_t i = 0; i < full_size*full_size; i++)
        copy_foc[i] = fock[i];

    dgeev_(&N,&V,&full_size,copy_foc,&full_size,eigRE,eigIM,eigL,&full_size,eigR,&full_size,work,&lwork,&info);
    // Sort MO's according to increasing energy
    sort_mos(full_size,eigRE,eigR);

    char T = 'T';
    double oned =1.;
    double zerod =0.;
 
    dgemm_(&N,&N,&full_size,&full_size,&full_size,&oned,sym_orth_mat,&full_size,eigR,&full_size,&zerod,c_mat,&full_size);

    //Memory leak while I get this working
    delete[] eigRE;
    delete[] eigIM;
    delete[] eigL;
    delete[] eigR;
    delete[] work;
    delete[] temp;
}

void hf_wfn::make_density_mat()
{
    int full_size = std::pow((-1+std::pow(1+8*mat_size,0.5))/2,1);
    density_mat = new double[full_size*full_size]{0.0};

    for(size_t i = 0; i < full_size; i++)
    {
        for(size_t j = 0; j < full_size; j++)
        {
            for(size_t m = 0; m < 5; m++)
            {
                density_mat[i*full_size + j] += c_mat[m*full_size+i]*c_mat[m*full_size+j];
            }
        }
    }
}

void hf_wfn::evaluate_energy()
{
    total_e=0.0;
    int nao = std::pow((-1+std::pow(1+8*mat_size,0.5))/2,1);
    double * ham_mat = new double[nao*nao]{0.0};
    double * result = new double[nao*nao]{0.0};
    
    char N = 'N';
    char T = 'T';
    double oned =1.;
    double zerod =0.;

    for(size_t i = 0; i < nao*nao; i++)
    {
        ham_mat[i]+=core_H[i]+inbasis_fock[i];
        //ham_mat[i]+=core_H[i]+result2[i];
    }

    dgemm_(&N,&N,&nao,&nao,&nao,&oned,density_mat,&nao,ham_mat,&nao,&zerod,result,&nao);

    for(size_t i = 0; i < nao; i++)
    {
        total_e+=result[i*nao+i];
    }
    std::cout << total_e << std::endl;
}

void hf_wfn::update_fock()
{
    int nao = std::pow((-1+std::pow(1+8*mat_size,0.5))/2,1);
    double * new_fock = new double[nao*nao]{0.0};
    double next_term;
    for(size_t i = 0; i < nao; i++)
    {
        for(size_t j = 0; j < nao; j++)
        {
            new_fock[i*nao + j] += core_H[i*nao + j];
            next_term=0.;
            for (size_t k = 0; k < nao; k++)
            {
                for(size_t l = 0; l < nao; l++)
                {
                    new_fock[i*nao + j] += density_mat[k*nao + l]*(2*teis(i+1,j+1,k+1,l+1)-teis(i+1,k+1,j+1,l+1));
                }
            }
        }
    }
    delete[]fock;

    fock=new_fock;
    for(size_t i = 0; i < nao*nao; i++)
        inbasis_fock[i]=fock[i];

    char N = 'N';
    char T = 'T';
    double oned =1.;
    double zerod =0.;

    double * temp = new double[nao*nao]{0.0};
    dgemm_(&T,&N,&nao,&nao,&nao,&oned,sym_orth_mat,&nao,fock,&nao,&zerod,temp,&nao);
    dgemm_(&N,&N,&nao,&nao,&nao,&oned,temp,&nao,sym_orth_mat,&nao,&zerod,fock,&nao);
    
}

void hf_wfn::print_mos()
{
    int nao = std::pow((-1+std::pow(1+8*mat_size,0.5))/2,1);
    print_mat(c_mat,nao);
}
