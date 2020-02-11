/* -*- indent-tabs-mode: t -*- */

#ifndef SOLVERS_POISSON
#define SOLVERS_POISSON

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

#include <math/complex.hpp>
#include <math/vec3d.hpp>
#include <math/array.hpp>
#include <multi/adaptors/fftw.hpp>
#include <basis/field.hpp>
#include <basis/fourier_space.hpp>

#include <fftw3-mpi.h>

namespace solvers {

	template<class basis_type>
  class poisson {

	public:

		auto operator()(const basis::field<basis_type, complex> & density) const {
			if(density.basis().periodic_dimensions() == 3){
				return solve_periodic(density);
			} else {
				return solve_finite(density);
			}
		}

		auto solve_periodic(const basis::field<basis_type, complex> & density) const {
			namespace fftw = boost::multi::fftw;
			
			basis::fourier_space fourier_basis(density.basis());

			basis::field<basis::fourier_space, complex> potential_fs(fourier_basis, density.basis_comm());

			if(not density.basis().dist().parallel()) {

				potential_fs.cubic() = fftw::dft(density.cubic(), fftw::forward);

			} else {
				
				ptrdiff_t local_n0, local_0_start;
				auto alloc_local = fftw_mpi_local_size_3d(density.basis().sizes()[0], density.basis().sizes()[1], density.basis().sizes()[2], &density.basis_comm(), &local_n0, &local_0_start);

				assert(local_n0 == density.basis().cubic_dist(0).local_size());
				assert(local_0_start == density.basis().cubic_dist(0).start());

				auto buff = fftw_alloc_complex(alloc_local);

				memcpy(buff, static_cast<const complex *>(density.linear().data()), density.linear().size()*sizeof(complex));
				
				auto plan = fftw_mpi_plan_dft_3d(density.basis().sizes()[0], density.basis().sizes()[1], density.basis().sizes()[2], buff, buff, &density.basis_comm(), FFTW_FORWARD, FFTW_ESTIMATE);

				fftw_execute(plan);

				memcpy(static_cast<complex *>(potential_fs.linear().data()), buff, potential_fs.linear().size()*sizeof(complex));

				fftw_free(buff);
							
				fftw_destroy_plan(plan);

			}

			const double scal = (-4.0*M_PI)/potential_fs.basis().size();
			
			for(int ix = 0; ix < potential_fs.basis().local_sizes()[0]; ix++){
				for(int iy = 0; iy < potential_fs.basis().local_sizes()[1]; iy++){
					for(int iz = 0; iz < potential_fs.basis().local_sizes()[2]; iz++){

						auto ixg = potential_fs.basis().cubic_dist(0).local_to_global(ix);
						auto iyg = potential_fs.basis().cubic_dist(1).local_to_global(iy);
						auto izg = potential_fs.basis().cubic_dist(2).local_to_global(iz);

						auto g2 = potential_fs.basis().g2(ixg, iyg, izg);

						if(potential_fs.basis().g_is_zero(ixg, iyg, izg) or fourier_basis.outside_sphere(g2)){
							potential_fs.cubic()[ix][iy][iz] = 0.0;
							continue;
						}
						
						potential_fs.cubic()[ix][iy][iz] *= -scal/g2;
					}
				}
			}

			basis::field<basis_type, complex> potential_rs(density.basis(), density.basis_comm());

			if(not density.basis().dist().parallel()) {

				potential_rs.cubic() = fftw::dft(potential_fs.cubic(), fftw::backward);
				
			} else {

				ptrdiff_t local_n0, local_0_start;

				auto alloc_local = fftw_mpi_local_size_3d(potential_fs.basis().sizes()[0], potential_fs.basis().sizes()[1], potential_fs.basis().sizes()[2], &potential_fs.basis_comm(), &local_n0, &local_0_start);
				
				assert(local_n0 == potential_fs.basis().cubic_dist(0).local_size());
				assert(local_0_start == potential_fs.basis().cubic_dist(0).start());

				auto buff = fftw_alloc_complex(alloc_local);

				memcpy(buff, static_cast<complex *>(potential_fs.linear().data()), potential_fs.linear().size()*sizeof(complex));
				
				auto plan = fftw_mpi_plan_dft_3d(potential_fs.basis().sizes()[0], potential_fs.basis().sizes()[1], potential_fs.basis().sizes()[2], buff, buff, &potential_fs.basis_comm(), FFTW_BACKWARD, FFTW_ESTIMATE);

				fftw_execute(plan);
				
				memcpy(static_cast<complex *>(potential_rs.linear().data()), buff, potential_rs.linear().size()*sizeof(complex));

				fftw_free(buff);
				
				fftw_destroy_plan(plan);

			}

			return potential_rs;
			
		}

		auto solve_finite(const basis::field<basis_type, complex> & density) const {
			namespace fftw = boost::multi::fftw;

			basis::field<basis_type, complex> potential2x(density.basis().enlarge(2));

			potential2x = 0.0;
			
			for(int ix = 0; ix < density.basis().sizes()[0]; ix++){
				for(int iy = 0; iy < density.basis().sizes()[1]; iy++){
					for(int iz = 0; iz < density.basis().sizes()[2]; iz++){
						auto i2x = ix;
						auto i2y = iy;
						auto i2z = iz;
						if(ix >= density.basis().sizes()[0]/2) i2x += density.basis().sizes()[0];
						if(iy >= density.basis().sizes()[1]/2) i2y += density.basis().sizes()[1];
						if(iz >= density.basis().sizes()[2]/2) i2z += density.basis().sizes()[2];
						potential2x.cubic()[i2x][i2y][i2z] = density.cubic()[ix][iy][iz];
					}
				}
			}
			
			fftw::dft_inplace(potential2x.cubic(), fftw::forward);

			basis::fourier_space fourier_basis(potential2x.basis());

			const auto scal = (-4.0*M_PI)/fourier_basis.size();
			const auto cutoff_radius = potential2x.basis().min_rlength()/2.0;

			for(int ix = 0; ix < fourier_basis.sizes()[0]; ix++){
				for(int iy = 0; iy < fourier_basis.sizes()[1]; iy++){
					for(int iz = 0; iz < fourier_basis.sizes()[2]; iz++){
						
						// this is the kernel of C. A. Rozzi et al., Phys. Rev. B 73, 205119 (2006).
						if(fourier_basis.g_is_zero(ix, iy, iz)){
							potential2x.cubic()[ix][iy][iz] *= -scal*cutoff_radius*cutoff_radius/2.0;
							continue;
						}
						
						auto g2 = fourier_basis.g2(ix, iy, iz);

						if(fourier_basis.outside_sphere(g2)){
							potential2x.cubic()[ix][iy][iz] = 0.0;
							continue;
						}

						potential2x.cubic()[ix][iy][iz] *= -scal*(1.0 - cos(cutoff_radius*sqrt(g2)))/g2;

					}
				}
			}

			fftw::dft_inplace(potential2x.cubic(), fftw::backward);
			
			basis::field<basis_type, complex> potential(density.basis());
			
			potential = 0.0;
			
			for(int ix = 0; ix < potential.basis().sizes()[0]; ix++){
				for(int iy = 0; iy < potential.basis().sizes()[1]; iy++){
					for(int iz = 0; iz < potential.basis().sizes()[2]; iz++){	
						auto i2x = ix;
						auto i2y = iy;
						auto i2z = iz;
						if(ix >= density.basis().sizes()[0]/2) i2x += density.basis().sizes()[0];
						if(iy >= density.basis().sizes()[1]/2) i2y += density.basis().sizes()[1];
						if(iz >= density.basis().sizes()[2]/2) i2z += density.basis().sizes()[2];
						potential.cubic()[ix][iy][iz] = potential2x.cubic()[i2x][i2y][i2z];
					}
				}
			}
			
			return potential;
		}
		
		basis::field<basis_type, double> operator()(const basis::field<basis_type, double> & density) const {

			using basis::field;
			
			//For the moment we copy to a complex array.
			field<basis_type, complex> complex_density(density.basis(), density.basis_comm());
			
			complex_density.linear() = density.linear();
			
			auto complex_potential = operator()(complex_density);
			
			field<basis_type, double> real_potential(density.basis(), density.basis_comm());
			
			//DATAOPERATIONS GPU::RUN 1D
			gpu::run(density.basis().dist().local_size(),
							 [rp = begin(real_potential.linear()), cp = begin(complex_potential.linear())]
							 GPU_LAMBDA (auto ii){
								 rp[ii] = real(cp[ii]);
							 });

			return real_potential;
			
		}
		
	private:
		
  };    
	
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>
#include <operations/integral.hpp>

TEST_CASE("class solvers::poisson", "[poisson]") {

	using namespace Catch::literals;
	namespace multi = boost::multi;
	using namespace basis;

	{

		fftw_mpi_init();

		auto comm = boost::mpi3::environment::get_world_instance();
		
		double ll = 10.0;
		
		ions::UnitCell cell({ll, 0.0, 0.0}, {0.0, ll, 0.0}, {0.0, 0.0, ll});
		basis::real_space rs(cell, input::basis::spacing(0.1), comm);

		SECTION("Grid periodic"){
		
			REQUIRE(rs.periodic_dimensions() == 3);
			
			REQUIRE(rs.sizes()[0] == 100);
			REQUIRE(rs.sizes()[1] == 100);
			REQUIRE(rs.sizes()[2] == 100);

		}
		
		field<real_space, complex> density(rs, comm);
		solvers::poisson<basis::real_space> psolver;
		
		SECTION("Point charge"){
		
			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						auto ixg = rs.cubic_dist(0).local_to_global(ix);
						auto iyg = rs.cubic_dist(1).local_to_global(iy);
						auto izg = rs.cubic_dist(2).local_to_global(iz);

						density.cubic()[ix][iy][iz] = 0.0;

						if(ixg == 0 and iyg == 0 and izg == 0) density.cubic()[ix][iy][iz] = -1.0;
					}
				}
			}
		
			auto potential = psolver(density);
		
			double sum[2] = {0.0, 0.0};
			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						sum[0] += fabs(real(potential.cubic()[ix][iy][iz]));
						sum[1] += fabs(imag(potential.cubic()[ix][iy][iz]));
					}
				}
			}

			comm.all_reduce_in_place_n(sum, 2, std::plus<>{});
			
			// These values haven't been validated against anything, they are
			// just for consistency. Of course the imaginary part has to be
			// zero, since the density is real.
		
			REQUIRE(sum[0] == 59.7758543176_a);
			REQUIRE(fabs(sum[1]) <= 1e-12);
		
			if(rs.cubic_dist(0).start() == 0 and rs.cubic_dist(1).start() == 0 and rs.cubic_dist(2).start() == 0) REQUIRE(real(potential.cubic()[0][0][0]) == -0.0241426581_a);
		}

		SECTION("Plane wave"){

			double kk = 2.0*M_PI/rs.rlength()[0];
		
			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						
						auto ixg = rs.cubic_dist(0).local_to_global(ix);
						auto iyg = rs.cubic_dist(1).local_to_global(iy);
						auto izg = rs.cubic_dist(2).local_to_global(iz);
						
						double xx = rs.rvector(ixg, iyg, izg)[0];
						density.cubic()[ix][iy][iz] = complex(cos(kk*xx), sin(kk*xx));
					
					}
				}
			}

			auto potential = psolver(density);

			double diff = 0.0;
			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						diff += fabs(potential.cubic()[ix][iy][iz] - 4*M_PI/kk/kk*density.cubic()[ix][iy][iz]);
					}
				}
			}

			comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});

			diff /= rs.size();
		
			REQUIRE(diff < 1.0e-14);
	
		}


		SECTION("Real plane wave"){

			field<real_space, double> rdensity(rs, comm);

			double kk = 8.0*M_PI/rs.rlength()[1];
		
			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						
						auto ixg = rs.cubic_dist(0).local_to_global(ix);
						auto iyg = rs.cubic_dist(1).local_to_global(iy);
						auto izg = rs.cubic_dist(2).local_to_global(iz);						
						double yy = rs.rvector(ixg, iyg, izg)[1];
						rdensity.cubic()[ix][iy][iz] = cos(kk*yy);
					}
				}
			}

			auto rpotential = psolver(rdensity);

			double diff = 0.0;
			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						diff += fabs(rpotential.cubic()[ix][iy][iz] - 4*M_PI/kk/kk*rdensity.cubic()[ix][iy][iz]);
					}
				}
			}

			diff /= rs.size();
		
			REQUIRE(diff < 1e-8);

		}
	}


	{

		const double ll = 8.0;
		
		ions::UnitCell cell({ll, 0.0, 0.0}, {0.0, ll, 0.0}, {0.0, 0.0, ll}, 0);
		basis::real_space rs(cell, input::basis::spacing(0.09));

		solvers::poisson<basis::real_space> psolver;

		SECTION("Grid finite"){		

			REQUIRE(rs.periodic_dimensions() == 0);
			
			REQUIRE(rs.sizes()[0] == 89);
			REQUIRE(rs.sizes()[1] == 89);
			REQUIRE(rs.sizes()[2] == 89);

		}
		
		field<real_space, complex> density(rs);

		SECTION("Point charge finite"){

			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						density.cubic()[ix][iy][iz] = 0.0;
						if(norm(rs.rvector(ix, iy, iz)) < 1e-10) density.cubic()[ix][iy][iz] = -1.0/rs.volume_element();
					}
				}
			}

			REQUIRE(real(operations::integral(density)) == -1.0_a);
			
			auto potential = psolver(density);

			std::ofstream ofile("pot.dat");
			
			//			double sumreal = 0.0;
			// double sumimag = 0.0;
			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						auto rr = rs.rvector(ix, iy, iz);
						if(iz == 0 and iy == 0) ofile << rr[0] << "\t" << real(potential.cubic()[ix][iy][iz]) << std::endl;						
					}
				}
			}

		}

	}

}


#endif


#endif
