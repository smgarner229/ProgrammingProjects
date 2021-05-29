#ifndef HF_WFN_HPP
#define HF_WFN_HPP

#include <map>
#include <utility>

class two_center_integral
{
    private:
        int i,j;
        double integral;
    public:
        double operator ()(int i, int j);
        bool operator ==(const two_center_integral & other);
};

class overlap_integral : public two_center_integral
{
    public:
        double operator()(int i, int j);
};

class kinetic_integral : public two_center_integral
{
    public:
        double operator()(int i, int j);
};



class hf_wfn
{
    double enuc;
};

#endif