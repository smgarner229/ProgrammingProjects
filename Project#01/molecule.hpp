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

    struct outofplane_angle
    {
        int i,j,k,l;
        double angle;

        outofplane_angle(int i, int j, int k, int l, double angle):
            i(i),j(j),k(k),l(l),angle(angle) {};
    };

    struct torsion_angle
    {
        int i,j,k,l;
        double angle;

        torsion_angle(int i, int j, int k, int l, double angle):
            i(i),j(j),k(k),l(l),angle(angle) {};
    };

    std::vector<particle> nuclei;
    std::vector<bond> bonds;
    std::vector<internal_angle> internal_angles;
    std::vector<outofplane_angle> outofplane_angles;
    std::vector<torsion_angle> torsion_angles;

    molecule(){}; //Constructor
    ~molecule(){};//Destructor

    void add_neucleus(double charge, double x, double y, double z);
    void calc_bond_lengths();
    void calc_bond_angles();
    void calc_outofplane_angle();
    

    friend std::ostream& operator <<(std::ostream & os, molecule & mol);

};

#endif