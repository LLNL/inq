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

#include <utils/profiling.hpp>

namespace inq {
namespace ions {

template <class cell_type, class geometry_type>
auto interaction_energy(const cell_type & cell, const geometry_type & geo, const hamiltonian::atomic_potential & atomic_pot){

	CALI_CXX_MARK_FUNCTION;

	double energy;
	boost::multi::array<math::vector3<double>, 1> forces(geo.num_atoms());
	boost::multi::array<double, 1> charges(geo.num_atoms());

	for(int ii = 0; ii < geo.num_atoms(); ii++) charges[ii] = atomic_pot.pseudo_for_element(geo.atoms()[ii]).valence_charge();

	interaction_energy(geo.num_atoms(), cell, charges, geo.coordinates(), atomic_pot.range_separation(), energy, forces);

	return energy;
}

template <class cell_type, class geometry_type>
auto interaction_forces(const cell_type & cell, const geometry_type & geo, const hamiltonian::atomic_potential & atomic_pot){

	CALI_CXX_MARK_FUNCTION;

	double energy;
	boost::multi::array<math::vector3<double>, 1> forces(geo.num_atoms());
	boost::multi::array<double, 1> charges(geo.num_atoms());

	for(int ii = 0; ii < geo.num_atoms(); ii++) charges[ii] = atomic_pot.pseudo_for_element(geo.atoms()[ii]).valence_charge();

	interaction_energy(geo.num_atoms(), cell, charges, geo.coordinates(), atomic_pot.range_separation(), energy, forces);

	return forces;
}

template <class cell_type, class array_charge, class array_positions, class array_forces>
void interaction_energy(const int natoms, const cell_type & cell, const array_charge & charge, const array_positions & positions, pseudo::math::erf_range_separation const & sep,
												double & energy, array_forces & forces){
	using math::vector3;

	const double alpha = 0.21;

	double ers = 0.0;
	for(int iatom = 0; iatom < natoms; iatom++) forces[iatom] = vector3<double>(0.0, 0.0, 0.0);

	double rcut = 6.0/alpha;

	for(int iatom = 0; iatom < natoms; iatom++){
		double zi = charge[iatom];

		double maxdist = 0.0;
		for(int jatom = 0; jatom < natoms; jatom++){
			maxdist = std::max(maxdist, norm(cell.position_in_cell(positions[jatom]) - cell.position_in_cell(positions[iatom])));
		}
		
		periodic_replicas rep(cell, cell.position_in_cell(positions[iatom]), rcut + sqrt(maxdist));

		for(unsigned irep = 0; irep < rep.size(); irep++){
			vector3<double> xi = rep[irep];
				
			for(int jatom = 0; jatom < natoms; jatom++){
				double zj = charge[jatom];
					
				vector3<double> rij = xi - cell.position_in_cell(positions[jatom]);
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
	for(int idir = 0; idir < 3; idir++) gcut = std::min(gcut, norm(cell.reciprocal(idir)));
	gcut = sqrt(gcut);
      
	const int isph = ceil(9.5*alpha/gcut);

	std::vector<std::complex<double> > phase(natoms);
    
	for(int ix = -isph; ix <= isph; ix++){
		for(int iy = -isph; iy <= isph; iy++){
			for(int iz = -isph; iz <= isph; iz++){
					
				const int ss = ix*ix + iy*iy + iz*iz;
					
				if(ss == 0 || ss > isph*isph) continue;
					
				vector3<double> gg = ix*cell.reciprocal(0) + iy*cell.reciprocal(1) + iz*cell.reciprocal(2);
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

#include <catch2/catch_all.hpp>

#include <vector>
#include <valarray>
#include <math/array.hpp>
#include <ions/unit_cell.hpp>

TEST_CASE("Function ions::interaction_energy", "[ions::interaction_energy]") {

	using namespace inq;
	using namespace Catch::literals;
  using math::vector3;
	const pseudo::math::erf_range_separation sep(0.625);
 
  SECTION("Aluminum cubic cell"){
  
    double aa = 7.653;
    
    ions::unit_cell cell(vector3<double>(aa, 0.0, 0.0), vector3<double>(0.0, aa, 0.0), vector3<double>(0.0, 0.0, aa));
    
    std::valarray<double> charge(4);
    charge = 3.0;
    
    std::vector<vector3<double>> positions(4);
    positions[0] = vector3<double>(0.0,    0.0,    0.0);
    positions[1] = vector3<double>(aa/2.0, aa/2.0, 0.0);
    positions[2] = vector3<double>(aa/2.0, 0.0,    aa/2.0);
    positions[3] = vector3<double>(0.0,    aa/2.0, aa/2.0);
    
    double energy;
    std::vector<vector3<double>> forces(4);
    
    ions::interaction_energy(4, cell, charge, positions, sep, energy, forces);
    
    CHECK(energy == -9.99517178_a); //this number comes from Octopus
    
  }

  SECTION("Diamond"){

    double aa = 6.74065308785213;
    
    ions::unit_cell cell(vector3<double>(0.0, aa/2.0, aa/2.0), vector3<double>(aa/2.0, 0.0, aa/2.0), vector3<double>(aa/2.0, aa/2.0, 0.0));

    const double charge[2] = {4.0, 4.0};
    
    std::vector<vector3<double>> positions(2);
    positions[0] = cell.metric().to_cartesian(vector3<double, math::contravariant>(0.0,  0.0,  0.0 ));
    positions[1] = cell.metric().to_cartesian(vector3<double, math::contravariant>(0.25, 0.25, 0.25));
    
    double energy;
    std::vector<vector3<double>> forces(2);

    ions::interaction_energy(2, cell, charge, positions, sep, energy, forces);

    CHECK(energy == -10.73490075_a); //this number comes from Octopus

  }

  SECTION("Iron"){
    
    double aa = 5.3970578;
    
    ions::unit_cell cell(vector3<double>(-aa/2.0, aa/2.0, aa/2.0), vector3<double>(aa/2.0, -aa/2.0, aa/2.0), vector3<double>(aa/2.0, aa/2.0, -aa/2.0));

    const double charge = 16.0;
    
    const vector3<double> position(0.0, 0.0, 0.0);
    
    double energy;
    std::vector<vector3<double>> forces(1);

    ions::interaction_energy(1, cell, &charge, &position, sep, energy, forces);

    CHECK(energy == -78.31680646_a); //this number comes from Octopus
    
  }
  
  SECTION("N2 supercell"){
    
    double aa = 20.0;
    
    ions::unit_cell cell(vector3<double>(aa, 0.0, 0.0), vector3<double>(0.0, aa, 0.0), vector3<double>(0.0, 0.0, aa));

		const double charge[2] = {5.0, 5.0};

    double distance = 2.0739744;

    std::vector<vector3<double>> positions(2);
    positions[0] = vector3<double>(0.0, 0.0, -distance/2.0);
    positions[1] = vector3<double>(0.0, 0.0,  distance/2.0);
    
    double energy;
    std::vector<vector3<double>> forces(2);

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

	
	SECTION("N2 supercell shifted"){
    
    double aa = 20.0;
    
    ions::unit_cell cell(vector3<double>(aa, 0.0, 0.0), vector3<double>(0.0, aa, 0.0), vector3<double>(0.0, 0.0, aa));

		const double charge[2] = {5.0, 5.0};

    double distance = 2.0739744;

    std::vector<vector3<double>> positions(2);
    positions[0] = vector3<double>(0.0, 0.0, 293.9 -distance/2.0);
    positions[1] = vector3<double>(0.0, 0.0, 293.9 + distance/2.0);
    
    double energy;
    std::vector<vector3<double>> forces(2);

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
		
	SECTION("Si64"){
    
    double alat = 20.25;
    
    ions::unit_cell cell(vector3<double>(alat, 0.0, 0.0), vector3<double>(0.0, alat, 0.0), vector3<double>(0.0, 0.0, alat));

		math::array<double, 1> charge(64, 4.0);

    std::vector<vector3<double>> positions;
		
		positions.push_back(math::vector3<double>(  alat*0.000000000,   alat*0.000000000,    alat*0.000000000));
	  positions.push_back(math::vector3<double>(  alat*0.000000000,   alat*0.000000000,    alat*0.500000000));
	  positions.push_back(math::vector3<double>(  alat*0.000000000,   alat*0.500000000,    alat*0.000000000));
	  positions.push_back(math::vector3<double>(  alat*0.000000000,   alat*0.500000000,    alat*0.500000000));
	  positions.push_back(math::vector3<double>(  alat*0.500000000,   alat*0.000000000,    alat*0.000000000));
	  positions.push_back(math::vector3<double>(  alat*0.500000000,   alat*0.000000000,    alat*0.500000000));
	  positions.push_back(math::vector3<double>(  alat*0.500000000,   alat*0.500000000,    alat*0.000000000));
	  positions.push_back(math::vector3<double>(  alat*0.500000000,   alat*0.500000000,    alat*0.500000000));
	  positions.push_back(math::vector3<double>(  alat*0.000000000,   alat*0.250000000,    alat*0.250000000));
	  positions.push_back(math::vector3<double>(  alat*0.000000000,   alat*0.250000000,    alat*0.750000000));
	  positions.push_back(math::vector3<double>(  alat*0.000000000,   alat*0.750000000,    alat*0.250000000));
	  positions.push_back(math::vector3<double>(  alat*0.000000000,   alat*0.750000000,    alat*0.750000000));
	  positions.push_back(math::vector3<double>(  alat*0.500000000,   alat*0.250000000,    alat*0.250000000));
	  positions.push_back(math::vector3<double>(  alat*0.500000000,   alat*0.250000000,    alat*0.750000000));
	  positions.push_back(math::vector3<double>(  alat*0.500000000,   alat*0.750000000,    alat*0.250000000));
	  positions.push_back(math::vector3<double>(  alat*0.500000000,   alat*0.750000000,    alat*0.750000000));
	  positions.push_back(math::vector3<double>(  alat*0.250000000,   alat*0.250000000,    alat*0.000000000));
	  positions.push_back(math::vector3<double>(  alat*0.250000000,   alat*0.250000000,    alat*0.500000000));
	  positions.push_back(math::vector3<double>(  alat*0.250000000,   alat*0.750000000,    alat*0.000000000));
	  positions.push_back(math::vector3<double>(  alat*0.250000000,   alat*0.750000000,    alat*0.500000000));
	  positions.push_back(math::vector3<double>(  alat*0.750000000,   alat*0.250000000,    alat*0.000000000));
	  positions.push_back(math::vector3<double>(  alat*0.750000000,   alat*0.250000000,    alat*0.500000000));
	  positions.push_back(math::vector3<double>(  alat*0.750000000,   alat*0.750000000,    alat*0.000000000));
	  positions.push_back(math::vector3<double>(  alat*0.750000000,   alat*0.750000000,    alat*0.500000000));
	  positions.push_back(math::vector3<double>(  alat*0.250000000,   alat*0.000000000,    alat*0.250000000));
	  positions.push_back(math::vector3<double>(  alat*0.250000000,   alat*0.000000000,    alat*0.750000000));
	  positions.push_back(math::vector3<double>(  alat*0.250000000,   alat*0.500000000,    alat*0.250000000));
	  positions.push_back(math::vector3<double>(  alat*0.250000000,   alat*0.500000000,    alat*0.750000000));
	  positions.push_back(math::vector3<double>(  alat*0.750000000,   alat*0.000000000,    alat*0.250000000));
	  positions.push_back(math::vector3<double>(  alat*0.750000000,   alat*0.000000000,    alat*0.750000000));
	  positions.push_back(math::vector3<double>(  alat*0.750000000,   alat*0.500000000,    alat*0.250000000));
	  positions.push_back(math::vector3<double>(  alat*0.750000000,   alat*0.500000000,    alat*0.750000000));
	  positions.push_back(math::vector3<double>(  alat*0.375000000,   alat*0.125000000,    alat*0.375000000));
	  positions.push_back(math::vector3<double>(  alat*0.375000000,   alat*0.125000000,    alat*0.875000000));
	  positions.push_back(math::vector3<double>(  alat*0.375000000,   alat*0.625000000,    alat*0.375000000));
	  positions.push_back(math::vector3<double>(  alat*0.375000000,   alat*0.625000000,    alat*0.875000000));
	  positions.push_back(math::vector3<double>(  alat*0.875000000,   alat*0.125000000,    alat*0.375000000));
	  positions.push_back(math::vector3<double>(  alat*0.875000000,   alat*0.125000000,    alat*0.875000000));
	  positions.push_back(math::vector3<double>(  alat*0.875000000,   alat*0.625000000,    alat*0.375000000));
	  positions.push_back(math::vector3<double>(  alat*0.875000000,   alat*0.625000000,    alat*0.875000000));
	  positions.push_back(math::vector3<double>(  alat*0.125000000,   alat*0.125000000,    alat*0.125000000));
	  positions.push_back(math::vector3<double>(  alat*0.125000000,   alat*0.125000000,    alat*0.625000000));
	  positions.push_back(math::vector3<double>(  alat*0.125000000,   alat*0.625000000,    alat*0.125000000));
	  positions.push_back(math::vector3<double>(  alat*0.125000000,   alat*0.625000000,    alat*0.625000000));
	  positions.push_back(math::vector3<double>(  alat*0.625000000,   alat*0.125000000,    alat*0.125000000));
	  positions.push_back(math::vector3<double>(  alat*0.625000000,   alat*0.125000000,    alat*0.625000000));
	  positions.push_back(math::vector3<double>(  alat*0.625000000,   alat*0.625000000,    alat*0.125000000));
	  positions.push_back(math::vector3<double>(  alat*0.625000000,   alat*0.625000000,    alat*0.625000000));
	  positions.push_back(math::vector3<double>(  alat*0.125000000,   alat*0.375000000,    alat*0.375000000));
	  positions.push_back(math::vector3<double>(  alat*0.125000000,   alat*0.375000000,    alat*0.875000000));
	  positions.push_back(math::vector3<double>(  alat*0.125000000,   alat*0.875000000,    alat*0.375000000));
	  positions.push_back(math::vector3<double>(  alat*0.125000000,   alat*0.875000000,    alat*0.875000000));
	  positions.push_back(math::vector3<double>(  alat*0.625000000,   alat*0.375000000,    alat*0.375000000));
	  positions.push_back(math::vector3<double>(  alat*0.625000000,   alat*0.375000000,    alat*0.875000000));
	  positions.push_back(math::vector3<double>(  alat*0.625000000,   alat*0.875000000,    alat*0.375000000));
	  positions.push_back(math::vector3<double>(  alat*0.625000000,   alat*0.875000000,    alat*0.875000000));
	  positions.push_back(math::vector3<double>(  alat*0.375000000,   alat*0.375000000,    alat*0.125000000));
	  positions.push_back(math::vector3<double>(  alat*0.375000000,   alat*0.375000000,    alat*0.625000000));
	  positions.push_back(math::vector3<double>(  alat*0.375000000,   alat*0.875000000,    alat*0.125000000));
	  positions.push_back(math::vector3<double>(  alat*0.375000000,   alat*0.875000000,    alat*0.625000000));
	  positions.push_back(math::vector3<double>(  alat*0.875000000,   alat*0.375000000,    alat*0.125000000));
	  positions.push_back(math::vector3<double>(  alat*0.875000000,   alat*0.375000000,    alat*0.625000000));
	  positions.push_back(math::vector3<double>(  alat*0.875000000,   alat*0.875000000,    alat*0.125000000));
	  positions.push_back(math::vector3<double>(  alat*0.875000000,   alat*0.875000000,    alat*0.625000000));
		
    double energy;
    std::vector<vector3<double>> forces(64);

    ions::interaction_energy(64, cell, charge, positions, sep, energy, forces);

    CHECK(energy == -253.0283966274_a); 

		for(int iatom = 0; iatom < 64; iatom++){
			CHECK(length(forces[iatom]) < 1e-12);
		}
  }
	
}
#endif

#endif
