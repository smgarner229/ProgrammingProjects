#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include <vector>
#include <iostream>

#define MAXSTRING 80

class particle
{
    public:

    double charge;
    double x,y,z;

    particle(){};
    particle(double charge, double x, double y, double z) : charge(charge),x(x),y(y),z(z) {std::cout << "added" << charge << " " << x << " " << y << " " << z<<std::endl;};
    ~particle(){};

    void print_ptcl()
    {
        std::cout << charge << " " << x << " " << y << " " << z << std::endl;
        return;
    };
};

class molecule
{
    public:
    
    std::vector<particle> nuclei;

    molecule(){}; //Constructor
    ~molecule(){};//Destructor

    void add_neucleus(double charge, double x, double y, double z)
    {   
        nuclei.push_back(particle(charge,x,y,z));
        return;
    }

    friend std::ostream& operator <<(std::ostream & os, molecule & mol)
    {
        os << "Charge\t\tx\t\ty\t\tz\n";
        for (int i = 0; i < mol.nuclei.size(); i++)
        {
            mol.nuclei[i].print_ptcl();
        }
        os << "\n";
    }

};

#endif