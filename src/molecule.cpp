#include "molecule.hpp"
#include "geometry_opps.hpp"

#define FIELDWIDTH 10

    // Friend class to print a point's coordinates
    std::ostream & operator <<(std::ostream & os, point & pt)
    {   
        os << std::fixed <<
              std::setw(FIELDWIDTH) << std::setprecision(5) << pt.x << 
              std::setw(FIELDWIDTH) << std::setprecision(5) << pt.y << 
              std::setw(FIELDWIDTH) << std::setprecision(5) << pt.z;
        return os;
    }

    // Friend class to print a particle.
    std::ostream & operator <<(std::ostream & os, particle & ptcl)
    {
        point temp(ptcl.x,ptcl.y,ptcl.z);
        os << std::setw(FIELDWIDTH) << std::setprecision(3) << ptcl.charge <<
              std::setw(FIELDWIDTH) << std::setprecision(4) << ptcl.mass   <<
              temp;
        return os;
    }

    // Translate a single particle by some distance
    void point::translate(double deltax, double deltay, double deltaz)
    {
        x+=deltax;
        y+=deltay;
        z+=deltaz;
    }

    // Constructor which initializes the inertial tensor to be all zeros
    molecule::molecule()
    {

    }

    // Destructor.  Clear memory allocated in the initializer
    molecule::~molecule()
    {
        delete inertial_tensor;
        delete hessian;
    }

    // Add a nucleus to the molecule.  
    // Mass is added via the particle constructor
    void molecule::add_neucleus(double charge, double x, double y, double z)
    {   
        nuclei.push_back(particle(charge,x,y,z));
        return;
    }

    // Friend class to print a molecule's geometry
    std::ostream& operator <<(std::ostream & os, molecule & mol)
    {
        os << std::setw(FIELDWIDTH) << "Charge" 
           << std::setw(FIELDWIDTH) << "Mass"  
           << std::setw(FIELDWIDTH) << "x"
           << std::setw(FIELDWIDTH) << "y"
           << std::setw(FIELDWIDTH) << "z" << std::endl;
        for (size_t i = 0; i < mol.nuclei.size(); i++)
            os << mol.nuclei[i] << std::endl;
        os << "\n";
        return os;
    }

    void molecule::calc_bond_lengths(bool print)
    {
        for(size_t i = 0; i < nuclei.size(); i++)
        {
            for(size_t j = i+1; j < nuclei.size(); j++)
            {
                bonds.push_back(bond(i,j,calc_distances(nuclei[i],nuclei[j])));
            }
        }
        if(print)
        {
            std::cout << "\nBond Distances:\n";
            for(size_t i = 0; i < bonds.size(); i++)
            {
                std::cout << bonds[i].to << " " << bonds[i].from << " " << bonds[i].dist << std::endl;
            }
        }
        return;
    }

    void molecule::calc_bond_angles(bool print)
    {
        for(size_t i = 0; i < nuclei.size(); i++)
        {
            for(size_t j = i+1; j < nuclei.size(); j++)
            {
                for(size_t k = j+1; k < nuclei.size(); k++)
                {
                    if(calc_distances(nuclei[i],nuclei[j])<4.0 && calc_distances(nuclei[j],nuclei[k])<4.0)
                    {
                        internal_angles.push_back(internal_angle(i,j,k,calc_angle(nuclei[i],nuclei[j],nuclei[k])));
                    }
                }
            }
        }
        if (print)
        {
            std::cout << "\nBond Angles:\n";
            for(size_t i=0; i < internal_angles.size(); i++)
            {
                std::cout << internal_angles[i].central << " " << internal_angles[i].to << 
                " " << internal_angles[i].from << " " << internal_angles[i].angle << std::endl;
            }
        }
    }

    void molecule::calc_outofplane_angle(bool print)
    {
        for(size_t i = 0; i < nuclei.size(); i++)
        {
            for(size_t j = 0; j < nuclei.size(); j++)
            {
                for(size_t k = 0; k < nuclei.size(); k++)
                {
                    for(size_t l = 0; l < j; l++)
                    {
                        if(calc_distances(nuclei[i],nuclei[k]) < 4.0 && calc_distances(nuclei[k],nuclei[j]) < 4.0 && calc_distances(nuclei[k],nuclei[l]) < 4.0
                            && i != j && j != k && k != l && i != k && i != l)
                        {
                            outofplane_angles.push_back(outofplane_angle(i,j,k,l,calc_out_of_plane_angle(nuclei[i],nuclei[j],nuclei[k],nuclei[l])));
                            if(print)
                            {
                                std::cout << i << " " << j << " " << k << " " << l << " " << outofplane_angles.at(outofplane_angles.size()-1).angle << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }


    void molecule::calc_torsion_angle(bool print)
    {
        for(size_t i = 0; i < nuclei.size(); i++)
        {
            for(size_t j = 0; j < i; j++)
            {
                for(size_t k = 0; k < j; k++)
                {
                    for(size_t l = 0; l < k; l++)
                    {
                        if(calc_distances(nuclei[i],nuclei[j]) < 4.0 && calc_distances(nuclei[j],nuclei[k]) < 4.0 && calc_distances(nuclei[k],nuclei[l]) < 4.0)
                        {
                            torsion_angles.push_back(torsion_angle(i,j,k,l,calc_torsional_angle(nuclei[i],nuclei[j],nuclei[k],nuclei[l])));
                            if(print)
                            {
                                std::cout << i << " " << j << " " << k << " " << l << " " << torsion_angles.at(torsion_angles.size()-1).angle << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }

    void molecule::calc_center_of_mass()
    {
        for(size_t i = 0; i < nuclei.size(); i++)
        {
            total_mass += nuclei[i].mass;
            com.x += nuclei[i].mass*nuclei[i].x;
            com.y += nuclei[i].mass*nuclei[i].y;
            com.z += nuclei[i].mass*nuclei[i].z;
        }
        com.x /= total_mass;
        com.y /= total_mass;
        com.z /= total_mass;
        std::cout << "\nCOM:\n";
        std::cout << com.x << "\t" << com.y << "\t" << com.z << std::endl;
        for(size_t i = 0; i < nuclei.size(); i++)
        {
            nuclei[i].translate(-com.x,-com.y,-com.z);
        }
    }

    void molecule::calc_inertial_tensor()
    {

        inertial_tensor = new double[9] {0.0};
        for(size_t i = 0; i < nuclei.size(); i++)
        {
            inertial_tensor[0] += nuclei[i].mass * (std::pow(nuclei[i].y,2.)+std::pow(nuclei[i].z,2.));
            inertial_tensor[4] += nuclei[i].mass * (std::pow(nuclei[i].x,2.)+std::pow(nuclei[i].z,2.));
            inertial_tensor[8] += nuclei[i].mass * (std::pow(nuclei[i].x,2.)+std::pow(nuclei[i].y,2.));

            inertial_tensor[1] += nuclei[i].mass * nuclei[i].x * nuclei[i].y;
            inertial_tensor[2] += nuclei[i].mass * nuclei[i].x * nuclei[i].z;
            inertial_tensor[5] += nuclei[i].mass * nuclei[i].y * nuclei[i].z;
        }
        inertial_tensor[3] = inertial_tensor[1];
        inertial_tensor[6] = inertial_tensor[2];
        inertial_tensor[7] = inertial_tensor[5];

        int n = 3;
        char Nchar='N';
        double *eigReal=new double[n];
        double *eigImag=new double[n];
        int one=1;
        int lwork=6*n;
        double *work=new double[lwork];
        int info;

        dgeev_(&Nchar,&Nchar,&n,inertial_tensor,&n,eigReal,eigImag,nullptr,&one,nullptr,&one,work,&lwork,&info); 

        std::cout << "Eigenvalues of Inertial Tensor:\n";
        for (size_t i = 0; i < n; i++)
        {
            std::cout << eigReal[i] << " + " << eigImag[i] << " i\n";
        }

        delete eigReal;
        delete eigImag;
        delete work;

    }

    void molecule::vibrational_analysis()
    {
        // Appropriately Mass Weight the read in Hessian
        for(size_t i = 0; i < 3*nuclei.size(); i++)
        {
            for(size_t j = 0; j < 3*nuclei.size(); j++)
            {
                hessian[i*3*nuclei.size()+j] = hessian[i*3*nuclei.size()+j] / std::pow(nuclei[i/3].mass * nuclei[j/3].mass,0.5);
            }
        }

        // Parameters to setup the LAPACK call
        int n = 3*nuclei.size();
        char Nchar='N';
        double *eigReal=new double[n];
        double *eigImag=new double[n];
        int one=1;
        int lwork=6*n;
        double *work=new double[lwork];
        int info;

        dgeev_(&Nchar,&Nchar,&n,hessian,&n,eigReal,eigImag,nullptr,&one,nullptr,&one,work,&lwork,&info); 

        std::sort(eigReal,eigReal+n);

        std::cout << "\n" << "Computed Eigen values of Hessian:\n";
        
        for (size_t i = 0; i < n; i++)
        {
            std::cout << eigReal[i] << " + " << eigImag[i] << " i\n";
        }

        delete eigReal;
        delete eigImag;
        delete work;
    }