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