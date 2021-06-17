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
    return _teis.find(_get_compound_index(i,j,k,l))->second;
}
