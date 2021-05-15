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
    particle(double charge, double x, double y, double z) : charge(charge),x(x),y(y),z(z) {};
    ~particle(){};

    void print_ptcl()
    {
        std::cout << charge << "\t\t" << x << "\t\t" << y << "\t\t" << z << std::endl;
        return;
    };
};

class molecule
{
    public:
    
    struct bond
    {
        public:
            int to;
            int from;
            double dist;

        bond(int to, int from, double dist):to(to),from(from),dist(dist){};
        ~bond(){};
    };

    struct internal_angle
    {
        public:
        int central;
        int to, from;
        double angle;

        internal_angle(int central, int to, int from, double angle):central(central), to(to), from(from), angle(angle){};
        ~internal_angle(){};
    };

    std::vector<particle> nuclei;
    std::vector<bond> bonds;
    std::vector<internal_angle> internal_angles;

    molecule(){}; //Constructor
    ~molecule(){};//Destructor

    void add_neucleus(double charge, double x, double y, double z);
    void calc_bond_lengths();
    void calc_bond_angles();

    friend std::ostream& operator <<(std::ostream & os, molecule & mol);

};

#endif