#ifndef TEI_HANDLER_HPP
#define TEI_HANDLER_HPP
#include <iostream>
#include <map>

class two_electron_integral_handler
{
    private:
        std::map<int,double> _teis;
        int _get_compound_index(int i, int j, int k, int l);
    public:
        two_electron_integral_handler(){};
        ~two_electron_integral_handler(){};
        
        // Ask for the tei element given the four indices we care about
        double operator()(int i, int j, int k, int l);

        // Add a piece to the tei's
        void add_tei(int i, int j, int k, int l, double tei);

        // Print all integrals
        void dump_tei();

        // Performs N^8 Orbital rotation to MO basis
        void slow_tei_transform(const double * c_mat, const int & norb);

        // Perform N^5 Orbital rotation to MO basis
        void rotate_integrals(const double * c_mat, const int & norb);

};


#endif