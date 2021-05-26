#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#ifndef PRINT
#define PRINT 0
#endif 

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>

#define NDIM 3

extern "C" {
    extern int dgeev_(char*,char*,int*,double*,int*,double*, double*, double*, int*, double*, int*, double*, int*, int*);
    extern int dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
}

const std::map <int,double> mass_map {{1,1.00797},{2,4.00260},{3,6.941},{4,9.01218},{5,10.81},{6,12.011},{7,14.0067},{8,15.9994},
                                      {9,18.9984403},{10,20.179}};

class point{
    public:
        double x,y,z;
    point(double x, double y, double z): x(x), y(y), z(z) {};
    point(){};
    ~point(){};

    void translate(double deltax, double deltay, double deltaz);
    friend std::ostream & operator << (std::ostream & os, point & pt);
};                              

class particle : public point
{
    public:

    int charge;
    
    double mass;

    particle(){};
    particle(int charge, double x, double y, double z) : point(x,y,z), charge(charge) {mass=mass_map.find(charge)->second;};
    ~particle(){};

    friend std::ostream & operator << (std::ostream & os, particle & ptcl);
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

    double total_mass;
    point com;
    double * inertial_tensor = nullptr;
    double * hessian = nullptr;

    molecule(); //Constructor
    ~molecule();//Destructor

    void add_neucleus(double charge, double x, double y, double z);
    void calc_bond_lengths(bool print = PRINT);
    void calc_bond_angles(bool print = PRINT);
    void calc_outofplane_angle(bool print = PRINT);
    void calc_torsion_angle(bool print = PRINT);

    void calc_center_of_mass();
    void calc_inertial_tensor();

    void vibrational_analysis();
    

    friend std::ostream& operator <<(std::ostream & os, molecule & mol);

};

#endif