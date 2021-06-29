#include "hf_wfn.hpp"
#include <cmath>
#include <iomanip>

#define FIELDWIDTH 10

extern "C" {
    extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
    extern int dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
    extern int dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
}

template <class T>
static void swap(T * one, T * two)
{
    T temp = *one;
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
            if(std::abs(mat[index])>1e-9) std::cout << std::fixed << std::setw(FIELDWIDTH) << std::setprecision(5) << mat[index];
            else std::cout << std::setw(FIELDWIDTH) << std::setprecision(5) << 0.0;
        }
    }
    std::cout << "\n\n";
    return;
}


double * triangle_to_full_mat(double * & triangle, const int & tri_size)
{
    double * full_mat = nullptr;
    int nao = std::pow((-1+std::pow(1+8*tri_size,0.5))/2,1);
    full_mat = new double[nao*nao];
    size_t counteri=0,counterj=0;

    for(size_t i = 0; i < tri_size; i++)
    {
        if (counteri == counterj)
        {
            full_mat[counteri*nao+counterj]=triangle[i];
            counteri++;
            counterj=0;
        }
        else
        {
            full_mat[counteri*nao+counterj]=triangle[i];
            full_mat[counterj*nao+counteri]=triangle[i];
            counterj++;
        }
    }
    return full_mat;
}

void hf_wfn::make_core_H()
{
    core_H = new double[mat_size];
    for(size_t j = 0; j < mat_size; j++)
    {
        core_H[j] = ke_ints[j] + eN_ints[j];
    }
    core_H = triangle_to_full_mat(core_H,mat_size);
    return;
}

void hf_wfn::orthogonalize_basis()
{
    double * temp = new double[nao*nao];
    temp = triangle_to_full_mat(sints,mat_size);

    char N = 'N';
    char V = 'V';
    char L = 'L';    
    char T = 'T';
    
    double * eigRE = new double[nao]{0.0};
    double * work = new double[4*nao]{0.0};
    double * orthmat = new double[nao*nao]{0.0};
    double * eigval_mat = new double[nao*nao]{0.0};
    sym_orth_mat = new double[nao*nao]{0.0};

    double oned =  1.;
    double zerod = 0.;
    int lwork = 4 * nao;
    int info;

    // Diagonalize symmetric overlap matrix
    // Note we're inputting the full matrix (even though dsyev only needs one of the triangles)
    // But it will be overwritten in the end anyway
    dsyev_(&V,&L,&nao,temp,&nao,eigRE,work,&lwork,&info);

    // Create the matrix with 1/sqrt(eigenvalues) on the diagonal
    for (size_t i = 0; i < nao; i++)
    { eigval_mat[i*nao+i] = std::pow(eigRE[i],-0.5); }

    // Reuse the temp matrix to store the result of L*S
    dgemm_(&N,&N,&nao,&nao,&nao,&oned,temp,&nao,eigval_mat,&nao,&zerod,orthmat,&nao);
    // Carry out (L*S)*L^T, store the result in sym_orth_mat
    dgemm_(&N,&T,&nao,&nao,&nao,&oned,orthmat,&nao,temp,&nao,&zerod,sym_orth_mat,&nao);

    delete[] temp;
    delete[] eigRE;
    delete[] work;
    delete[] eigval_mat;
    delete[] orthmat;
    return;
}

void hf_wfn::make_initial_fock()
{
    fock= new double[nao*nao]{0.0};
    double * temp = new double[nao*nao]{0.0};

    char N = 'N';
    char T = 'T';
    double oned =1.;
    double zerod =0.;

    inbasis_fock = new double[nao*nao]{0.0};
    for(size_t i = 0; i < nao*nao; i++)
        inbasis_fock[i]=core_H[i];

    dgemm_(&T,&N,&nao,&nao,&nao,&oned,sym_orth_mat,&nao,core_H,&nao,&zerod,temp,&nao);
    dgemm_(&N,&N,&nao,&nao,&nao,&oned,temp,&nao,sym_orth_mat,&nao,&zerod,fock,&nao);

    delete[] temp;
    return;
}

static void sort_mos(int nmo, double * & mo_energy, double * & mos)
{
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
    double * temp = new double[nao*nao];
    double * eigRE = new double[nao]{0.0};
    double * work = new double[nao*nao]{0.0};
    c_mat = new double[nao*nao];
    
    char N = 'N';
    char V = 'V';
    char L = 'L';

    int lwork = 6 * nao;
    int info = 0;

    double oned =1.;
    double zerod =0.;    
    
    // Diagonalize fock matrix (which is symmetric), and sort the resulting eigenvalues & eigenvectors
    // \mathbf{F}_0^'\mathbf{C}_0^'=\mathbf{C}_0^'\mathbf{\epsilon}_0
    dsyev_(&V,&L,&nao,fock,&nao,eigRE,work,&lwork,&info);
    sort_mos(nao,eigRE,fock);
 
    if(!orb_energies) delete[] orb_energies;
    orb_energies = new double[nao];
    for(size_t i = 0; i < nao; i++)
        orb_energies[i] = eigRE[i];

    // \mathbf{S}^{-1/2}\mathbf{C_0^'}
    dgemm_(&N,&N,&nao,&nao,&nao,&oned,sym_orth_mat,&nao,fock,&nao,&zerod,c_mat,&nao);

    delete[] eigRE;
    delete[] work;

    return;
}

void hf_wfn::make_density_mat()
{
    if(niter) // Not the first iteration
    {
        delete[] last_density_mat;
        last_density_mat = new double[nao*nao];
        for(size_t i = 0; i < nao*nao; i++)
        {
            last_density_mat[i] = density_mat[i];
        }
        delete[] density_mat;
    }
    else // First iteration
    {
        last_density_mat = new double[nao*nao]{0.0};
        niter++;
    }
    density_mat = new double[nao*nao]{0.0};

    for(size_t i = 0; i < nao; i++)
    {
        for(size_t j = 0; j < nao; j++)
        {
            for(size_t m = 0; m < 5; m++)
            {
                density_mat[i*nao + j] += c_mat[m*nao+i]*c_mat[m*nao+j];
            }
        }
    }
}

void hf_wfn::evaluate_energy()
{
    last_e  = total_e;
    total_e = enuc;
    double * ham_mat = new double[nao*nao]{0.0};
    double * result = new double[nao*nao]{0.0};
    
    char N = 'N';
    char T = 'T';
    double oned =1.;
    double zerod =0.;

    for(size_t i = 0; i < nao*nao; i++)
    {ham_mat[i]+=core_H[i]+inbasis_fock[i];}

    dgemm_(&N,&N,&nao,&nao,&nao,&oned,density_mat,&nao,ham_mat,&nao,&zerod,result,&nao);

    for(size_t i = 0; i < nao; i++)
    {total_e+=result[i*nao+i];}

    delete[] ham_mat;
    delete[] result;
    return;
}

void hf_wfn::update_fock()
{
    // Hold onto the new and old fock matrix, to see if we've converged
    double * new_fock = new double[nao*nao]{0.0};
    // For matrix multiplication
    double * temp = new double[nao*nao]{0.0};

    for(size_t i = 0; i < nao; i++)
    {
        for(size_t j = 0; j < nao; j++)
        {
            new_fock[i*nao + j] += core_H[i*nao + j];
            for (size_t k = 0; k < nao; k++)
            {
                for(size_t l = 0; l < nao; l++)
                {
                    new_fock[i*nao + j] += density_mat[k*nao + l]*(2*teis(i+1,j+1,k+1,l+1)-teis(i+1,k+1,j+1,l+1));
                }
            }
        }
    }
    // Compare new fock matrix to old fock matrix?
    delete[]fock;

    fock=new_fock;
    for(size_t i = 0; i < nao*nao; i++)
        inbasis_fock[i]=fock[i];

    char N = 'N';
    char T = 'T';
    double oned =1.;
    double zerod =0.;

    dgemm_(&T,&N,&nao,&nao,&nao,&oned,sym_orth_mat,&nao,fock,&nao,&zerod,temp,&nao);
    dgemm_(&N,&N,&nao,&nao,&nao,&oned,temp,&nao,sym_orth_mat,&nao,&zerod,fock,&nao);
    
    delete[] temp;
    return;
}

void hf_wfn::print_mos()
{
    print_mat(c_mat,nao,true);
}

void hf_wfn::fock_in_mo_basis()
{
    double * new_fock = new double[nao*nao]{0.0};
    for(size_t i = 0; i < nao; i++)
    {
        for(size_t j = 0; j < nao; j++)
        {
            for(size_t k = 0; k < nao; k++)
            {
                for(size_t l = 0; l < nao; l++)
                {
                    new_fock[i*nao+j]+=c_mat[i*nao+k]*c_mat[j*nao+l]*inbasis_fock[i*nao + j];
                }
            }
        }
    }
    print_mat(new_fock,nao);
}

bool hf_wfn::check_converged()
{
    double dmdiff = 0.0;
    for(size_t i = 0; i < nao*nao; i++)
    dmdiff+=std::pow(density_mat[i]-last_density_mat[i],2);
    dmdiff = std::pow(dmdiff,0.5);
    std::cout << std::fixed << std::setprecision(15) << total_e <<  "\t" << std::abs(last_e-total_e) << "\t" << dmdiff << std::endl;
    return dmdiff < 0.00000000001 && std::abs(last_e-total_e) < 0.00000000001;
}