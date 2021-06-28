#include <iostream>

#include "tei_handler.hpp"

static inline void swap(int & a, int & b)
{
    a = a + b;
    b = a - b;
    a = a - b;
    return;
}

int two_electron_integral_handler::_get_compound_index(int i, int j, int k, int l)
{
    if (j > i) swap(i,j);
    if (l > k) swap(l,k);
    int ij = i*(i+1)/2+j, kl = k*(k+1)/2+l;
    if (kl>ij) swap(ij,kl);
    return ij*(ij+1)/2+kl;
}

void two_electron_integral_handler::add_tei(int i, int j, int k, int l, double tei)
{
    _teis.insert(std::pair<int,double>(_get_compound_index(i,j,k,l),tei));
}

void two_electron_integral_handler::dump_tei()
{
    for (std::map<int,double>::iterator i = _teis.begin(); i != _teis.end(); i++)
    {
        std::cout << i->first << "\t" << i->second << std::endl;
    }
    return;
}

double two_electron_integral_handler::operator()(int i, int j, int k, int l)
{
    std::map<int,double>::iterator eri_val = _teis.find(_get_compound_index(i,j,k,l));
    if (eri_val != _teis.end())
        return eri_val->second;
    return 0.0;
}

void two_electron_integral_handler::slow_tei_transform(const double * c_mat, const int & norb)
{
    std::map<int, double> holder;
    for(size_t i = 0; i < norb; i++)
    {
        for(size_t j = 0; j <= i; j++)
        {
            for(size_t k = 0; k <= i; k++)
            {
                for(size_t l = 0; l <= (i==k ? j : k); l++)
                {
                    double new_integral = 0.0;
                    for(size_t p = 0; p < norb; p++)
                    {
                        for(size_t q = 0; q < norb; q++)
                        {
                            for(size_t r = 0; r < norb; r++)
                            {
                                for(size_t s = 0; s < norb; s++)
                                {
                                    new_integral+=(*this)(p+1,q+1,r+1,s+1)*c_mat[i*norb+p]*c_mat[j*norb+q]*c_mat[k*norb+r]*c_mat[l*norb+s];
                                }
                            }
                        }
                    }
                    holder.insert(std::pair<int,double>(_get_compound_index(i+1,j+1,k+1,l+1),new_integral));
                }
            }
        }
    }
    this->_teis=holder;
}

void two_electron_integral_handler::rotate_integrals(const double * c_mat, const int & norb)
{
    std::map<int, double> holder;
    for(size_t i = 0; i < norb; i++)
    {
        double integral = 0.0;
        for(size_t p = 0; p < norb; p++)
        {
            for(size_t q = 0; q < norb; q++)
            {
                for(size_t r = 0; r < norb; r++)
                {
                    for(size_t s = 0; s < norb; s++)
                    {
                        integral += c_mat[i*norb + p] * (*this)(p+1,q+r,r+1,s+1);
                    }
                }
            }
        }
    }
}