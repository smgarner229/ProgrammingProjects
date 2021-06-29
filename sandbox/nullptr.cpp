#include <iostream>

int main(int argc, char ** argv)
{

    double * dp = nullptr;
    double * tp = nullptr;
    tp = new double(0.0);
    if(dp)
    {
        std::cout << "Nullptr did print\n";
    }
    else
    {
        std::cout << "Nullptr DID NOT print\n";
    }

    if(tp)
    {
        std::cout << "Nullptr did print\n";
    }
    else
    {
        std::cout << "Nullptr DID NOT print\n";
    }



}
