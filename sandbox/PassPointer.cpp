#include <iostream>
#include <algorithm>

double * send_new_ptr(size_t nds);

int main(int argc, char ** argv)
{
    double * holder = nullptr;
    holder = send_new_ptr(25);
    for(size_t i = 0; i < 25; i++)
        std::cout << i << "\t" << holder[i] << std::endl;
    delete[] holder;
    holder = send_new_ptr(10);
    delete[] holder;
    for(size_t i = 0; i < 25; i++)
        std::cout << i << "\t" << holder[i] << std::endl;
    return 0;
}

double * send_new_ptr(size_t nds)
{
    double * return_value = nullptr;
    return_value = new double[nds]{1.0};
    std::fill_n(return_value,nds,-1);
    return return_value;
}