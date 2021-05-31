#include <iostream>

int main(int argc, char** argv)
{
    double * one = nullptr;
    double * two = nullptr;
    one = new double[9]{0.0};
    one[3]=3.;
    two = one;
    std::cout << two[3];
    delete[] one;
    std::cout << two[3];
}