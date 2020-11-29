/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__IONS__INTERACTION
#define INQ__IONS__INTERACTION

/*
 Copyright (C) 2019 Xavier Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <inq_config.h>
#include <math/vector3.hpp>
#include <cmath>
#include <ions/periodic_replicas.hpp>
#include <multi/array.hpp>
#include <limits>
#include <complex>
#include <hamiltonian/atomic_potential.hpp>

#include <caliper/cali.h>

namespace inq {
namespace ions {

template <class cell_type, class geometry_type>
auto interaction_energy(const cell_type & cell, const geometry_type & geo, const hamiltonian::atomic_potential & atomic_pot){

	CALI_CXX_MARK_FUNCTION;

	double energy;
	boost::multi::array<math::vec3d, 1> forces(geo.num_atoms());
	boost::multi::array<double, 1> charges(geo.num_atoms());

	for(int ii = 0; ii < geo.num_atoms(); ii++) charges[ii] = atomic_pot.pseudo_for_element(geo.atoms()[ii]).valence_charge();

	interaction_energy(geo.num_atoms(), cell, charges, geo.coordinates(), atomic_pot.range_separation(), energy, forces);

	return energy;
}

template <class cell_type, class geometry_type>
auto interaction_forces(const cell_type & cell, const geometry_type & geo, const hamiltonian::atomic_potential & atomic_pot){

	CALI_CXX_MARK_FUNCTION;

	double energy;
	boost::multi::array<math::vec3d, 1> forces(geo.num_atoms());
	boost::multi::array<double, 1> charges(geo.num_atoms());

	for(int ii = 0; ii < geo.num_atoms(); ii++) charges[ii] = atomic_pot.pseudo_for_element(geo.atoms()[ii]).valence_charge();

	interaction_energy(geo.num_atoms(), cell, charges, geo.coordinates(), atomic_pot.range_separation(), energy, forces);

	return forces;
}

template <class cell_type, class array_charge, class array_positions, class array_forces>
void interaction_energy(const int natoms, const cell_type & cell, const array_charge & charge, const array_positions & positions, pseudo::math::erf_range_separation const & sep,
												double & energy, array_forces & forces){
	using math::vec3d;

	const double alpha = 0.21;
		
	double ers = 0.0;
	for(int iatom = 0; iatom < natoms; iatom++) forces[iatom] = vec3d(0.0, 0.0, 0.0);

	double rcut = 6.0/alpha;

	for(int iatom = 0; iatom < natoms; iatom++){
		double zi = charge[iatom];
      
		periodic_replicas rep(cell, positions[iatom], rcut);

		for(unsigned irep = 0; irep < rep.size(); irep++){
			vec3d xi = rep[irep];
				
			for(int jatom = 0; jatom < natoms; jatom++){
				double zj = charge[jatom];
					
				vec3d rij = xi - positions[jatom];
				double rr = length(rij);
					
				if(rr < 1.0e-5) continue;
					
				auto eor = erfc(alpha*rr)/rr;
					
				ers += 0.5*zi*zj*eor;

				forces[jatom] -= zi*zj*rij*(eor + 2.0*alpha/sqrt(M_PI)*exp(-alpha*alpha*rr*rr))/(rr*rr);
			}
		}
      
	}

	double eself = 0.0;

	// self-interaction
	double total_charge = 0.0;
	for(int iatom = 0; iatom < natoms; iatom++){
		double zi = charge[iatom];

		total_charge += zi;
		eself -= alpha*zi*zi/sqrt(M_PI);
	}

	// G = 0 energy
	auto efs = -M_PI*total_charge*total_charge/(2.0*alpha*alpha*cell.volume());

	double gcut = std::numeric_limits<double>::max();
	for(int idir = 0; idir < 3; idir++) gcut = std::min(gcut, norm(cell.b(idir)));
	gcut = sqrt(gcut);
      
	const int isph = ceil(9.5*alpha/gcut);

	std::vector<std::complex<double> > phase(natoms);
    
	for(int ix = -isph; ix <= isph; ix++){
		for(int iy = -isph; iy <= isph; iy++){
			for(int iz = -isph; iz <= isph; iz++){
					
				const int ss = ix*ix + iy*iy + iz*iz;
					
				if(ss == 0 || ss > isph*isph) continue;
					
				vec3d gg = ix*cell.b(0) + iy*cell.b(1) + iz*cell.b(2);
				double gg2 = norm(gg);
					
				double exparg = -0.25*gg2/(alpha*alpha);
					
				if(exparg < -36.0) continue;
					
				double factor = 2.0*M_PI/cell.volume()*exp(exparg)/gg2;
					
				std::complex<double> sumatoms = 0.0;
				for(int iatom = 0; iatom < natoms; iatom++){
					double gx = dot(gg, positions[iatom]);
					auto aa = charge[iatom]*std::complex<double>(cos(gx), sin(gx));
					phase[iatom] = aa;
					sumatoms += aa;
				}
					
				efs += factor*std::real(sumatoms*std::conj(sumatoms));
					
				for(int iatom = 0; iatom < natoms; iatom++){
					for(int idir = 0; idir < 3; idir++){
						auto tmp = std::complex<double>(0.0, 1.0)*gg[idir]*phase[iatom];
						forces[iatom][idir] -= factor*std::real(std::conj(tmp)*sumatoms + tmp*std::conj(sumatoms));
					}
				}
					
			}
		}
	}

	// Previously unaccounted G = 0 term from pseudopotentials. 
	// See J. Ihm, A. Zunger, M.L. Cohen, J. Phys. C 12, 4409 (1979)
	double epseudo = 0.0;
	for(int iatom = 0; iatom < natoms; iatom++){
		epseudo += M_PI*charge[iatom]*pow(sep.sigma()*sqrt(2.0), 2)/cell.volume()*total_charge;
	}

	energy = ers + eself + efs + epseudo;

}
}
}

#ifdef INQ_IONS_INTERACTION_UNIT_TEST
#undef INQ_IONS_INTERACTION_UNIT_TEST

#include <catch2/catch.hpp>

#include <vector>
#include <valarray>
#include <ions/unitcell.hpp>

TEST_CASE("Function ions::interaction_energy", "[ions::interaction_energy]") {

	using namespace inq;
	using namespace Catch::literals;
  using math::vec3d;
	const pseudo::math::erf_range_separation sep(0.625);
 
  SECTION("Aluminum cubic cell"){
  
    double aa = 7.653;
    
    ions::UnitCell cell(vec3d(aa, 0.0, 0.0), vec3d(0.0, aa, 0.0), vec3d(0.0, 0.0, aa));
    
    std::valarray<double> charge(4);
    charge = 3.0;
    
    std::vector<vec3d> positions(4);
    positions[0] = vec3d(0.0,    0.0,    0.0);
    positions[1] = vec3d(aa/2.0, aa/2.0, 0.0);
    positions[2] = vec3d(aa/2.0, 0.0,    aa/2.0);
    positions[3] = vec3d(0.0,    aa/2.0, aa/2.0);
    
    double energy;
    std::vector<vec3d> forces(4);
    
    ions::interaction_energy(4, cell, charge, positions, sep, energy, forces);
    
    CHECK(energy == -9.99517178_a); //this number comes from Octopus
    
  }

  SECTION("Diamond"){

    double aa = 6.74065308785213;
    
    ions::UnitCell cell(vec3d(0.0, aa/2.0, aa/2.0), vec3d(aa/2.0, 0.0, aa/2.0), vec3d(aa/2.0, aa/2.0, 0.0));

    const double charge[2] = {4.0, 4.0};
    
    std::vector<vec3d> positions(2);
    positions[0] = cell.crystal_to_cart(vec3d(0.0,  0.0,  0.0 ));
    positions[1] = cell.crystal_to_cart(vec3d(0.25, 0.25, 0.25));
    
    double energy;
    std::vector<vec3d> forces(2);

    ions::interaction_energy(2, cell, charge, positions, sep, energy, forces);

    CHECK(energy == -10.73490075_a); //this number comes from Octopus

  }

  SECTION("Iron"){
    
    double aa = 5.3970578;
    
    ions::UnitCell cell(vec3d(-aa/2.0, aa/2.0, aa/2.0), vec3d(aa/2.0, -aa/2.0, aa/2.0), vec3d(aa/2.0, aa/2.0, -aa/2.0));

    const double charge = 16.0;
    
    const vec3d position(0.0, 0.0, 0.0);
    
    double energy;
    std::vector<vec3d> forces(1);

    ions::interaction_energy(1, cell, &charge, &position, sep, energy, forces);

    CHECK(energy == -78.31680646_a); //this number comes from Octopus
    
  }
  
  SECTION("N2 supercell"){
    
    double aa = 20.0;
    
    ions::UnitCell cell(vec3d(aa, 0.0, 0.0), vec3d(0.0, aa, 0.0), vec3d(0.0, 0.0, aa));

		const double charge[2] = {5.0, 5.0};

    double distance = 2.0739744;

    std::vector<vec3d> positions(2);
    positions[0] = vec3d(0.0, 0.0, -distance/2.0);
    positions[1] = vec3d(0.0, 0.0,  distance/2.0);
    
    double energy;
    std::vector<vec3d> forces(2);

    ions::interaction_energy(2, cell, charge, positions, sep, energy, forces);

		//these numbers come from Octopus
    CHECK(energy == 5.02018926_a); 

		CHECK(fabs(forces[0][0]) < 1.0e-16);
		CHECK(fabs(forces[0][1]) < 1.0e-16);
		CHECK(forces[0][2] == -5.7840844208_a);
		CHECK(fabs(forces[1][0]) < 1.0e-16);
		CHECK(fabs(forces[1][1]) < 1.0e-16);
		CHECK(forces[1][2] == 5.7840844208_a);
    
  }
	
}
#endif

#endif
