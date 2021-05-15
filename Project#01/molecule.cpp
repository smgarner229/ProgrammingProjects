#include "molecule.hpp"
#include "geometry_opps.hpp"

    void molecule::add_neucleus(double charge, double x, double y, double z)
    {   
        nuclei.push_back(particle(charge,x,y,z));
        return;
    }

    std::ostream& operator <<(std::ostream & os, molecule & mol)
    {
        os << "Charge\t\tx\t\ty\t\tz\n";
        for (size_t i = 0; i < mol.nuclei.size(); i++)
        {
            mol.nuclei[i].print_ptcl();
        }
        os << "\n";
        return os;
    }

    void molecule::calc_bond_lengths()
    {
        for(size_t i = 0; i < nuclei.size(); i++)
        {
            for(size_t j = i+1; j < nuclei.size(); j++)
            {
                bonds.push_back(bond(i,j,calc_distances(nuclei[i],nuclei[j])));
            }
        }
        for(size_t i = 0; i < bonds.size(); i++)
        {
            std::cout << bonds[i].to << " " << bonds[i].from << " " << bonds[i].dist << std::endl;
        }
        return;
    }

    void molecule::calc_bond_angles()
    {
        for(size_t i = 0; i < nuclei.size(); i++)
        {
            for(size_t j = i+1; j < nuclei.size(); j++)
            {
                for(size_t k = j+1; k < nuclei.size(); k++)
                {
                    if(calc_distances(nuclei[i],nuclei[j])<4.0 && calc_distances(nuclei[j],nuclei[k])<4.0)
                    {
                        internal_angles.push_back(internal_angle(j,i,k,calc_angle(nuclei[j],nuclei[i],nuclei[k])));
                    }
                }
            }
        }
        for(size_t i=0; i < internal_angles.size(); i++)
        {
            std::cout << internal_angles[i].central << " " << internal_angles[i].to << 
            " " << internal_angles[i].from << " " << internal_angles[i].angle << std::endl;
        }
    }