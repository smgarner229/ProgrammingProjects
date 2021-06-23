#include <iostream>
#include <bitset>

int main(int argc, char**argv)
{
    short int i=11,j=22,k=33,l=44;

    unsigned long ijkl=0;

    ijkl |= i;
    ijkl <<= 16;
    ijkl |= j;
    ijkl <<= 16;
    ijkl |= k;
    ijkl <<= 16;
    ijkl |= l;

    std::bitset<16> bi(i),bj(j),bk(k),bl(l);

    std::cout << bi << "\t" << bj << "\t" << bk << "\t" << bl << std::endl;
    std::bitset<64> z(ijkl);
    std::cout << z << "\t" << ijkl <<std::endl;

    return 0;
}
