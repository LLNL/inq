/* -*- indent-tabs-mode: t -*- */

#ifndef HAMILTONIAN_PROJECTOR
#define HAMILTONIAN_PROJECTOR

#include <multi/array.hpp>

#include <math/d3vector.hpp>
#include <math/spherical_harmonic.hpp>
#include <pseudo/pseudopotential.hpp>
#include <ions/unitcell.hpp>
#include <ions/periodic_replicas.hpp>
#include <basis/real_space.hpp>
#include <basis/spherical_grid.hpp>
#include <multi/adaptors/blas.hpp>

namespace hamiltonian {

  class projector {

  public:
    projector(const basis::real_space & basis, const ions::UnitCell & cell, pseudo::pseudopotential ps, math::d3vector atom_position):
      sphere_(basis, cell, atom_position, ps.projector_radius()),
      nproj_(ps.num_projectors_lm()),
      matrix_({nproj_, sphere_.size()}){

			std::vector<double> grid(sphere_.size()), proj(sphere_.size());

			int iproj_lm = 0;
      for(int iproj_l = 0; iproj_l < ps.num_projectors_l(); iproj_l++){
				// interpolate the value of the radial part of the projectors to the sphere points
				ps.projector(iproj_l).value(sphere_.size(), sphere_.distance(), proj);
				
				int l = ps.projector_l(iproj_l);
				
				// now construct the projector with the spherical harmonics
				for(int m = -l; m <= l; m++){
					for(int ipoint = 0; ipoint < sphere_.size(); ipoint++){
						matrix_[iproj_lm][ipoint] = proj[ipoint]*math::spherical_harmonic(l, m, sphere_.point_pos()[ipoint]);
					}
					iproj_lm++;
					kb_coeff_.push_back(ps.kb_coeff(iproj_l));
				}
				
      }

			assert(iproj_lm == ps.num_projectors_lm());
			
    }

    template <class field_set_type>
    void operator()(const field_set_type & phi, field_set_type & vnlphi) const {
			
      boost::multi::array<states::ks_states::coeff_type, 2> sphere_phi({sphere_.size(), phi.set_size()});

      sphere_.gather(phi.cubic(), sphere_phi);

      boost::multi::array<states::ks_states::coeff_type, 2> projections({nproj_, phi.set_size()});

			//DATAOPERATIONS BLAS
			boost::multi::blas::gemm('N', 'N', sphere_.volume_element(), sphere_phi, matrix_, 0.0, projections);
			
			//DATAOPERATIONS LOOP 2D
      for(int iproj = 0; iproj < nproj_; iproj++) {
				for(int ist = 0; ist < phi.set_size(); ist++)	projections[iproj][ist] *= kb_coeff_[iproj];
			}

			//DATAOPERATIONS BLAS
			boost::multi::blas::gemm('N', 'T', 1.0, projections, matrix_, 0.0, sphere_phi);
			
      sphere_.scatter_add(sphere_phi, vnlphi.cubic());
			
    }

    int num_projectors() const {
      return nproj_;
    }
		
    auto kb_coeff(int iproj){
      return kb_coeff_[iproj];
    }
		
  private:

    basis::spherical_grid sphere_;
    int nproj_;
		//OPTIMIZATION: make this matrix real
    boost::multi::array<complex, 2> matrix_;
    std::vector<double> kb_coeff_;
    
  };
  
}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class hamiltonian::projector", "[hamiltonian::projector]") {
  
  using namespace Catch::literals;
  using math::d3vector;
	
	const math::erf_range_separation sep(0.625);
	
	pseudo::pseudopotential ps(config::path::unit_tests_data() + "N.upf", sep);

  double ecut = 20.0;
  double ll = 10.0;

	ions::geometry geo;
  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  basis::real_space rs(cell, input::basis::cutoff_energy(ecut));

	hamiltonian::projector proj(rs, cell, ps, d3vector(0.0, 0.0, 0.0));

	REQUIRE(proj.num_projectors() == 8);
	
	REQUIRE(proj.kb_coeff(0) ==  7.494508815_a);
	REQUIRE(proj.kb_coeff(1) ==  0.6363049519_a);
	REQUIRE(proj.kb_coeff(2) == -4.2939052122_a);
	REQUIRE(proj.kb_coeff(3) == -4.2939052122_a);
	REQUIRE(proj.kb_coeff(4) == -4.2939052122_a);
	REQUIRE(proj.kb_coeff(5) == -1.0069878791_a);
	REQUIRE(proj.kb_coeff(6) == -1.0069878791_a);
	REQUIRE(proj.kb_coeff(7) == -1.0069878791_a);
	
}
#endif

#endif

