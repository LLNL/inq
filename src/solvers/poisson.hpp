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
#include <math/vector3.hpp>
#include <math/array.hpp>
#include <basis/field.hpp>
#include <basis/fourier_space.hpp>
#include <operations/space.hpp>
#include <operations/transfer.hpp>

#include <utils/profiling.hpp>

namespace inq {
namespace solvers {

class poisson {

public:

	basis::field<basis::real_space, complex> poisson_solve_periodic(const basis::field<basis::real_space, complex> & density) const {

		CALI_CXX_MARK_FUNCTION;
		
		const basis::real_space & real_space = density.basis();
		basis::fourier_space fourier_basis(real_space);

		auto potential_fs = operations::space::to_fourier(density);
			
		const double scal = (-4.0*M_PI)/fourier_basis.size();

		{
			CALI_CXX_MARK_SCOPE("poisson_finite_kernel_periodic");
			
			//DATAOPERATIONS GPU::RUN 4D
			gpu::run(fourier_basis.local_sizes()[2], fourier_basis.local_sizes()[1], fourier_basis.local_sizes()[0],
							 [point_op = fourier_basis.point_op(), pfs = begin(potential_fs.cubic()), scal] GPU_LAMBDA (auto iz, auto iy, auto ix){
								 
								 auto g2 = point_op.g2(ix, iy, iz);
								 
								 if(point_op.g_is_zero(ix, iy, iz)){
									 pfs[ix][iy][iz] = complex(0.0, 0.0);
								 } else {
									 pfs[ix][iy][iz] = pfs[ix][iy][iz]*(-scal/g2);
								 }
								 
							 });
		}
		
		return operations::space::to_real(potential_fs,  /*normalize = */ false);
	}

	basis::field<basis::real_space, complex> poisson_solve_finite(const basis::field<basis::real_space, complex> & density) const {

		CALI_CXX_MARK_FUNCTION;

		auto potential2x = operations::transfer::enlarge(density, density.basis().enlarge(2));
		auto potential_fs = operations::space::to_fourier(potential2x);
			
		auto fourier_basis = potential_fs.basis();

		const auto scal = (-4.0*M_PI)/fourier_basis.size();
		const auto cutoff_radius = potential2x.basis().min_rlength()/2.0;

		{
			CALI_CXX_MARK_SCOPE("poisson_finite_kernel_finite");

			gpu::run(fourier_basis.local_sizes()[2], fourier_basis.local_sizes()[1], fourier_basis.local_sizes()[0],
							 [point_op = fourier_basis.point_op(), pfs = begin(potential_fs.cubic()), scal, cutoff_radius] GPU_LAMBDA (auto iz, auto iy, auto ix){
								 
								 // this is the kernel of C. A. Rozzi et al., Phys. Rev. B 73, 205119 (2006).
								 if(point_op.g_is_zero(ix, iy, iz)){
									 pfs[ix][iy][iz] = pfs[ix][iy][iz]*(-scal)*cutoff_radius*cutoff_radius/2.0;
									 return;
								 }
								 
								 auto g2 = point_op.g2(ix, iy, iz);
								 pfs[ix][iy][iz] = pfs[ix][iy][iz]*(-scal)*(1.0 - cos(cutoff_radius*sqrt(g2)))/g2;
							 });
		}

		potential2x = operations::space::to_real(potential_fs,  /*normalize = */ false);
		auto potential = operations::transfer::shrink(potential2x, density.basis());

		return potential;
	}

	auto operator()(const basis::field<basis::real_space, complex> & density) const {

		CALI_CXX_MARK_SCOPE("poisson(complex)");
		
		if(density.basis().periodic_dimensions() == 3){
			return poisson_solve_periodic(density);
		} else {
			return poisson_solve_finite(density);
		}
	}
		
	basis::field<basis::real_space, double> operator()(const basis::field<basis::real_space, double> & density) const {

		CALI_CXX_MARK_SCOPE("poisson(real)");
		
		auto complex_potential = operator()(complex_field(density));
		return real_field(complex_potential);
	}

private:
		
};    
	
}
}


#ifdef INQ_SOLVERS_POISSON_UNIT_TEST
#undef INQ_SOLVERS_POISSON_UNIT_TEST

#include <catch2/catch.hpp>
#include <basis/real_space.hpp>
#include <ions/unitcell.hpp>
#include <operations/integral.hpp>

TEST_CASE("class solvers::poisson", "[solvers::poisson]") {


	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	namespace multi = boost::multi;
	using namespace basis;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	{

		systems::box box = systems::box::orthorhombic(10.0_b, 10.0_b, 13.7_b).spacing(0.1_b);
		basis::real_space rs(box, comm);

		SECTION("Grid periodic"){
		
			CHECK(rs.periodic_dimensions() == 3);
			
			CHECK(rs.sizes()[0] == 100);
			CHECK(rs.sizes()[1] == 100);
			CHECK(rs.sizes()[2] == 137);

		}
		
		field<real_space, complex> density(rs);
		solvers::poisson psolver;
		
		SECTION("Point charge"){
		
			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						auto ixg = rs.cubic_dist(0).local_to_global(ix);
						auto iyg = rs.cubic_dist(1).local_to_global(iy);
						auto izg = rs.cubic_dist(2).local_to_global(iz);

						density.cubic()[ix][iy][iz] = 0.0;

						if(ixg.value() == 0 and iyg.value() == 0 and izg.value() == 0) density.cubic()[ix][iy][iz] = -1.0;
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
		
			CHECK(sum[0] == 82.9383793318_a);
			CHECK(fabs(sum[1]) <= 5e-12);
		
			if(rs.cubic_dist(0).start() == 0 and rs.cubic_dist(1).start() == 0 and rs.cubic_dist(2).start() == 0) CHECK(real(potential.cubic()[0][0][0]) == -0.0241804443_a);
		}

		SECTION("Plane wave"){

			double kk = 2.0*M_PI/rs.rlength()[0];
		
			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						
						auto ixg = rs.cubic_dist(0).local_to_global(ix);
						auto iyg = rs.cubic_dist(1).local_to_global(iy);
						auto izg = rs.cubic_dist(2).local_to_global(iz);
						
						double xx = rs.point_op().rvector(ixg, iyg, izg)[0];
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
		
			CHECK(diff < 1.0e-13);
	
		}


		SECTION("Real plane wave"){

			field<real_space, double> rdensity(rs);

			double kk = 8.0*M_PI/rs.rlength()[1];
		
			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						
						auto ixg = rs.cubic_dist(0).local_to_global(ix);
						auto iyg = rs.cubic_dist(1).local_to_global(iy);
						auto izg = rs.cubic_dist(2).local_to_global(iz);						
						double yy = rs.point_op().rvector(ixg, iyg, izg)[1];
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

			comm.all_reduce_in_place_n(&diff, 1, std::plus<>{});

			diff /= rs.size();
		
			CHECK(diff < 1e-8);

		}
	}


	{

		systems::box box = systems::box::cubic(8.0_b).spacing(0.09_b).finite();
		basis::real_space rs(box, comm);

		solvers::poisson psolver;

		SECTION("Grid finite"){		

			CHECK(rs.periodic_dimensions() == 0);
			
			CHECK(rs.sizes()[0] == 89);
			CHECK(rs.sizes()[1] == 89);
			CHECK(rs.sizes()[2] == 89);

		}
		
		field<real_space, complex> density(rs);

		SECTION("Point charge finite"){

			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						density.cubic()[ix][iy][iz] = 0.0;
						if(norm(rs.point_op().rvector(ix, iy, iz)) < 1e-10) density.cubic()[ix][iy][iz] = -1.0/rs.volume_element();
					}
				}
			}

			CHECK(real(operations::integral(density)) == -1.0_a);
			
			auto potential = psolver(density);

			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
					for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
						
						auto ixg = rs.cubic_dist(0).local_to_global(ix);
						auto iyg = rs.cubic_dist(1).local_to_global(iy);
						auto izg = rs.cubic_dist(2).local_to_global(iz);

						auto rr = rs.point_op().rvector(ixg, iyg, izg).length();

						// it should be close to -1/r
						if(rr > 1) CHECK(fabs(potential.cubic()[ix][iy][iz]*rr + 1.0) < 0.025);
					
					}
				}
			}
			
			auto & part = potential.basis().part();
			if(part.contains(0))      CHECK(real(potential.linear()[part.global_to_local(utils::global_index(0))])      == -27.175214167_a);
			if(part.contains(8102))   CHECK(real(potential.linear()[part.global_to_local(utils::global_index(8102))])   ==  -2.9731998189_a);
			if(part.contains(700102)) CHECK(real(potential.linear()[part.global_to_local(utils::global_index(700102))]) ==  -0.2524115517_a);
			if(part.contains(27848))  CHECK(real(potential.linear()[part.global_to_local(utils::global_index(27848))])  ==  -0.2470080223_a);
			if(part.contains(612909)) CHECK(real(potential.linear()[part.global_to_local(utils::global_index(612909))]) ==  -0.2275712710_a);
			if(part.contains(368648)) CHECK(real(potential.linear()[part.global_to_local(utils::global_index(368648))]) ==  -0.1844298173_a);
			
			/*
			std::ofstream ofile("pot.dat");
			
			//			double sumreal = 0.0;
			// double sumimag = 0.0;
			for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
				auto rr = rs.rvector(ix, 0, 0);
				ofile << rr[0] << "\t" << real(potential.cubic()[ix][0][0]) << std::endl;						
			}
			*/

		}

	}

}


#endif


#endif
