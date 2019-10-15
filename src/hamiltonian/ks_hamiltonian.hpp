/* -*- indent-tabs-mode: t -*- */

#ifndef HAMILTONIAN_KS_HAMILTONIAN
#define HAMILTONIAN_KS_HAMILTONIAN

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

#include <basis/field.hpp>
#include <states/ks_states.hpp>
#include <multi/adaptors/fftw.hpp>
#include <hamiltonian/projector.hpp>
#include <operations/space.hpp>

namespace hamiltonian {
  template <class basis_type>
  class ks_hamiltonian {
		
  public:

		basis::field<basis::real_space, double> scalar_potential;
		boost::multi::array<double, 1> hf_occupations;
		basis::field_set<basis::real_space, complex> hf_orbitals;
		
		template <class orbitals_type>
    ks_hamiltonian(const basis_type & basis, const ions::UnitCell & cell, const atomic_potential & pot, const ions::geometry & geo,
									 const orbitals_type & orbitals):
			scalar_potential(basis),
			hf_occupations(orbitals.set_size()),
			hf_orbitals(basis, orbitals.set_size()) {

			scalar_potential = pot.local_potential(basis, cell, geo);

			for(int iatom = 0; iatom < geo.num_atoms(); iatom++){
				projectors_.push_back(projector(basis, cell, pot.pseudo_for_element(geo.atoms()[iatom]), geo.coordinates()[iatom]));
			}

    }

		void non_local(const basis::field_set<basis::real_space, complex> & phi, basis::field_set<basis::real_space, complex> & vnlphi) const {
			for(unsigned iproj = 0; iproj < projectors_.size(); iproj++) projectors_[iproj](phi, vnlphi);
		}

		auto non_local(const basis::field_set<basis::real_space, complex> & phi) const {
			basis::field_set<basis::real_space, complex> vnlphi(phi.basis(), phi.set_size());
			vnlphi = 0.0;
			non_local(phi, vnlphi);
			return vnlphi;
		}

		void fourier_space_terms(const basis::field_set<basis::fourier_space, complex> & phi, basis::field_set<basis::fourier_space, complex> & hphi) const {
			
			for(int ix = 0; ix < hphi.basis().gsize()[0]; ix++){
				for(int iy = 0; iy < hphi.basis().gsize()[1]; iy++){
					for(int iz = 0; iz < hphi.basis().gsize()[2]; iz++){
						double lapl = -0.5*(-hphi.basis().g2(ix, iy, iz));
						for(int ist = 0; ist < hphi.set_size(); ist++) hphi.cubic()[ix][iy][iz][ist] += lapl*phi.cubic()[ix][iy][iz][ist];
					}
				}
			}
			
		}
		
		void fourier_space_terms(basis::field_set<basis::fourier_space, complex> & hphi) const {
			
			for(int ix = 0; ix < hphi.basis().gsize()[0]; ix++){
				for(int iy = 0; iy < hphi.basis().gsize()[1]; iy++){
					for(int iz = 0; iz < hphi.basis().gsize()[2]; iz++){
						double lapl = -0.5*(-hphi.basis().g2(ix, iy, iz));
						for(int ist = 0; ist < hphi.set_size(); ist++) hphi.cubic()[ix][iy][iz][ist] *= lapl;
					}
				}
			}
			
		}
		
		void real_space_terms(const basis::field_set<basis::real_space, complex> & phi, basis::field_set<basis::real_space, complex> & hphi) const {
			//the non local potential in real space
			non_local(phi, hphi);

			//the scalar local potential in real space
			for(long ip = 0; ip < phi.basis().size(); ip++){
				double vv  = scalar_potential[ip];
				for(int ist = 0; ist < phi.set_size(); ist++) hphi[ip][ist] += vv*phi[ip][ist];
			}

			// Hartree-Fock exchange
			for(int ii = 0; ii < phi.set_size(); ii++){
				for(int jj = 0; jj < phi.set_size(); jj++){

					basis::field<basis::real_space, complex> rhoij(phi.basis());

					//DATAOPERATIONS LOOP 1D
					for(long ipoint = 0; ipoint < phi.basis().size(); ipoint++) rhoij[ipoint] = hf_orbitals[ipoint][jj]*phi[ipoint][ii];

					poisson_solver_(rhoij);

					//DATAOPERATIONS LOOP 1D
					for(long ipoint = 0; ipoint < phi.basis().size(); ipoint++) hphi[ipoint][ii] += hf_occupations[jj]*hf_orbitals[ipoint][jj]*rhoij[ipoint];
					
				}
			}
			
		}

    auto operator()(const basis::field_set<basis::real_space, complex> & phi) const{
      
			auto hphi_fs = operations::space::to_fourier(phi);

			fourier_space_terms(hphi_fs);
	
			auto hphi = operations::space::to_real(hphi_fs);

			real_space_terms(phi, hphi);

			return hphi;
			
		}
		
    auto operator()(const basis::field_set<basis::fourier_space, complex> & phi) const{

			auto phi_rs = operations::space::to_real(phi);

			basis::field_set<basis::real_space, complex> hphi_rs(phi_rs.basis(), phi.set_size());

			hphi_rs = 0.0;
						
			real_space_terms(phi_rs, hphi_rs);
		
			auto hphi = operations::space::to_fourier(hphi_rs);

			fourier_space_terms(phi, hphi);

			return hphi;
			
		}
		
		int num_projectors() const {
			int nn = 0;
			for(unsigned iproj = 0; iproj < projectors_.size(); iproj++){
				nn += projectors_[iproj].num_projectors();
			}
			return nn;			
		}

    template <class output_stream>
    void info(output_stream & out) const {
      out << "HAMILTONIAN:" << std::endl;
      out << "  Total number of projectors = " << num_projectors() << std::endl;
      out << std::endl;
    }	
		
  private:

		solvers::poisson<basis::real_space> poisson_solver_;
		std::vector<projector> projectors_;

  };

}

#ifdef UNIT_TEST

#include <ions/unitcell.hpp>
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>

TEST_CASE("Class hamiltonian::ks_hamiltonian", "[hamiltonian::ks_hamiltonian]"){

  using namespace Catch::literals;
  using math::d3vector;
  
  double ecut = 20.0;
  double ll = 10.0;

	ions::geometry geo;
  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  basis::real_space rs(cell, input::basis::cutoff_energy(ecut));

	SECTION("Basis"){
		
		REQUIRE(rs.size() == 8000);
		REQUIRE(rs.rspacing()[0] == 0.5_a);
		REQUIRE(rs.rspacing()[1] == 0.5_a);	
		REQUIRE(rs.rspacing()[2] == 0.5_a);
		REQUIRE(rs.volume_element() == 0.125_a);
	}
	
	hamiltonian::atomic_potential pot(geo.num_atoms(), geo.atoms());
	
	states::ks_states st(states::ks_states::spin_config::UNPOLARIZED, 11.0);

  basis::field_set<basis::real_space, complex> phi(rs, st.num_states());
	basis::field_set<basis::real_space, complex> hphi(rs, st.num_states());
	
	hamiltonian::ks_hamiltonian<basis::real_space> ham(rs, cell, pot, geo, phi);

	SECTION("Constant function"){
		
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] = 1.0;
					}
				}
			}
		}
		
		hphi = ham(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++){
						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - 0.0);
					}
				}
			}
		}

		diff /= hphi.cubic().num_elements();
		
		REQUIRE(diff == 0.0_a);
		
	}
	
	SECTION("Plane wave"){
		
		double kk = 2.0*M_PI/rs.rlength()[0];
		
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++){
						double xx = rs.rvector(ix, iy, iz)[0];
						phi.cubic()[ix][iy][iz][ist] = complex(cos(ist*kk*xx), sin(ist*kk*xx));
					}
				}
			}
		}

		hphi = ham(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++){
							diff += fabs(hphi.cubic()[ix][iy][iz][ist] - 0.5*ist*kk*ist*kk*phi.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		diff /= hphi.cubic().num_elements();

		REQUIRE(diff < 1e-15);
		
	}


	SECTION("Harmonic oscillator"){

		double ww = 2.0;

		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					double r2 = rs.r2(ix, iy, iz);
					ham.scalar_potential.cubic()[ix][iy][iz] = ww*r2;

					for(int ist = 0; ist < phi.set_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] = exp(-0.5*ww*r2);
					}
					
				}
			}
		}

		hphi = ham(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++){
						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - ww*phi.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		diff /= hphi.cubic().num_elements();

		REQUIRE(diff == 0.0055687279_a);
		
	}


	SECTION("Plane wave - fourier"){
		
		double kk = 2.0*M_PI/rs.rlength()[0];
		
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++){
						double xx = rs.rvector(ix, iy, iz)[0];
						phi.cubic()[ix][iy][iz][ist] = complex(cos(ist*kk*xx), sin(ist*kk*xx));
					}
				}
			}
		}

		hphi = operations::space::to_real(ham(operations::space::to_fourier(phi)));
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++){
							diff += fabs(hphi.cubic()[ix][iy][iz][ist] - 0.5*ist*kk*ist*kk*phi.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		diff /= hphi.cubic().num_elements();

		REQUIRE(diff < 1e-15);
		
	}


	SECTION("Harmonic oscillator"){

		double ww = 2.0;

		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					double r2 = rs.r2(ix, iy, iz);
					ham.scalar_potential.cubic()[ix][iy][iz] = ww*r2;

					for(int ist = 0; ist < phi.set_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] = exp(-0.5*ww*r2);
					}
					
				}
			}
		}

		hphi = ham(phi);
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++){
						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - ww*phi.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		diff /= hphi.cubic().num_elements();

		REQUIRE(diff == 0.0055687279_a);
		
	}


	SECTION("Harmonic oscillator - fourier"){

		double ww = 2.0;

		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					double r2 = rs.r2(ix, iy, iz);
					ham.scalar_potential.cubic()[ix][iy][iz] = ww*r2;

					for(int ist = 0; ist < phi.set_size(); ist++){
						phi.cubic()[ix][iy][iz][ist] = exp(-0.5*ww*r2);
					}
					
				}
			}
		}

		hphi = operations::space::to_real(ham(operations::space::to_fourier(phi)));
		
		double diff = 0.0;
		for(int ix = 0; ix < rs.rsize()[0]; ix++){
			for(int iy = 0; iy < rs.rsize()[1]; iy++){
				for(int iz = 0; iz < rs.rsize()[2]; iz++){
					for(int ist = 0; ist < phi.set_size(); ist++){
						diff += fabs(hphi.cubic()[ix][iy][iz][ist] - ww*phi.cubic()[ix][iy][iz][ist]);
					}
				}
			}
		}

		diff /= hphi.cubic().num_elements();

		REQUIRE(diff == 0.0055687279_a);
		
	}
	
}

#endif

#endif
