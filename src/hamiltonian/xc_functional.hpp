/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__XC_FUNCTIONAL
#define INQ__HAMILTONIAN__XC_FUNCTIONAL

/*
 Copyright (C) 2019-2023 Xavier Andrade, Alexey Kartsev

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

#include <stdexcept>

#include <xc.h>
#include <observables/density.hpp>
#include <operations/gradient.hpp>
#include <operations/integral.hpp>
#include <operations/divergence.hpp>
#include <operations/laplacian.hpp>
#include <basis/field.hpp>

#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace hamiltonian {
	class xc_functional {

		public:
		
		auto true_functional() const {
			return id_ > 0;			
		}
		
		xc_functional(const int functional_id, int spin_components):
			id_(functional_id),
			nspin_(spin_components){

			assert(nspin_ == 1 or nspin_ == 2);
			
			if(true_functional() and xc_func_init(&func_, functional_id, (nspin_==1)?XC_UNPOLARIZED:XC_POLARIZED) != 0){
				fprintf(stderr, "Functional '%d' not found\n", functional_id);
				exit(1);
			}
		}

		xc_functional(xc_functional const & other):
			xc_functional(other.id_, other.nspin_){
		}
		
		xc_functional operator=(xc_functional const & other) {
			if(id_ == other.id_ and nspin_ == other.nspin_) return *this;
			if(true_functional()) xc_func_end(&func_);
			id_ = other.id_;
			nspin_ = other.nspin_;
			if(true_functional() and xc_func_init(&func_, id_, (nspin_==1)?XC_UNPOLARIZED:XC_POLARIZED) != 0){
				fprintf(stderr, "Functional '%d' not found\n", id_);
				exit(1);
			}
			return *this;
		}
		
		~xc_functional(){
			if(true_functional()) xc_func_end(&func_);
		}

		template <class field_type>
		void operator()(field_type const & spin_density, double & xc_energy, field_type & vxc) const {

			assert(true_functional());
			
			CALI_CXX_MARK_SCOPE("xc_functional");
			
			basis::field<basis::real_space, double>  exc(vxc.basis());
			
			assert(spin_density.matrix().num_elements() == vxc.matrix().num_elements());
			
			switch(func_.info->family) {
				case XC_FAMILY_LDA:{
					xc_lda_exc_vxc(&func_, spin_density.basis().local_size(), raw_pointer_cast(spin_density.matrix().data_elements()), raw_pointer_cast(exc.linear().data_elements()), raw_pointer_cast(vxc.matrix().data_elements()));
					gpu::sync();
					break;
				}
				case XC_FAMILY_GGA:{
					ggafunctional(spin_density, exc, vxc);
					break;
				}	
				default:{
					std::runtime_error("inq error: unsupported xc functional family");
					break;
				}
			}

			xc_energy = operations::integral_product(exc, observables::density::total(spin_density));
		}

		auto & libxc_func() const {
			return func_;
		}
		
		auto libxc_func_ptr() const {
			return &func_;
		}

		auto exx_coefficient() const {
			if(xc_hyb_type(&func_) == XC_HYB_HYBRID) return xc_hyb_exx_coef(&func_);
			return 0.0;
		}
		
		//this function has to be public because of cuda limitations
		template <class density_type, class exc_type, class vxc_type>
		void ggafunctional(density_type const & density, exc_type & exc, vxc_type & vxc) const { 

			// Info about this implementation:
			//
			//   https://tddft.org/programs/libxc/manual/libxc-5.1.x/
			//   
			
			CALI_CXX_MARK_FUNCTION;

			auto grad_real = operations::gradient(density);

			auto nsigma = 1;
			if(nspin_ > 1) nsigma = 3;
			
			basis::field_set<basis::real_space, double> sigma(vxc.basis(), nsigma);

			gpu::run(vxc.basis().local_size(),
							 [sig = begin(sigma.matrix()), grad = begin(grad_real.matrix()), metric = density.basis().cell().metric(), nsigma] GPU_LAMBDA (auto ip){
								 sig[ip][0] = metric.norm(grad[ip][0]);
								 if(nsigma > 1) {
									 sig[ip][1] = metric.dot(grad[ip][0], grad[ip][1]);
									 sig[ip][2] = metric.norm(grad[ip][1]);
								 }
							 });

			basis::field_set<basis::real_space, double> vsigma(vxc.basis(), nsigma);
			
			xc_gga_exc_vxc(&func_, density.basis().local_size(), density.data(), sigma.data(), exc.data(), vxc.data(), vsigma.data());
			gpu::sync();
			
			basis::field_set<basis::real_space, math::vector3<double, math::covariant>> vxc_extra(vxc.skeleton());

			gpu::run(vxc.basis().local_size(),
							 [vex = begin(vxc_extra.matrix()), vsig = begin(vsigma.matrix()), grad = begin(grad_real.matrix()), nsigma] GPU_LAMBDA (auto ip){
								 if(nsigma == 1) {
										vex[ip][0] = 2.0*vsig[ip][0]*grad[ip][0];
								 } else {
										vex[ip][0] = 2.0*vsig[ip][0]*grad[ip][0] + vsig[ip][1]*grad[ip][1];
										vex[ip][1] = 2.0*vsig[ip][2]*grad[ip][1] + vsig[ip][1]*grad[ip][0];
								 }
							 });

			auto divvxcextra = operations::divergence(vxc_extra);

			gpu::run(vxc.local_set_size(), vxc.basis().local_size(),
							 [vv = begin(vxc.matrix()), div = begin(divvxcextra.matrix())] GPU_LAMBDA (auto ispin, auto ip){
								 vv[ip][ispin] -= div[ip][ispin];
							 });
		}

	private:
		
		int id_;
		int nspin_;
		xc_func_type func_;
			
	};

}
}


///////////////////////////////////////////////////////////////////

#ifdef INQ_HAMILTONIAN_XC_FUNCTIONAL_UNIT_TEST
#undef INQ_HAMILTONIAN_XC_FUNCTIONAL_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <operations/randomize.hpp>
#include <operations/gradient.hpp>
#include <basis/field.hpp>
#include <math/array.hpp>
#include <utils/finite_difference.hpp>
#include <ions/geometry.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>

GPU_FUNCTION auto sqwave(inq::math::vector3<double> rr, int n){
	auto kvec = 2.0*M_PI*inq::math::vector3<double>(1.0/9.0, 1.0/12.0, 1.0/10.0);
	return  sin(n*dot(kvec, rr))*sin(n*dot(kvec, rr));
}

auto laplacian_sqwave(inq::math::vector3<double> rr, int n){
	auto kvec = 2.0*M_PI*inq::math::vector3<double>(1.0/9.0, 1.0/12.0, 1.0/10.0);
	auto ff = 0.0;
	for(int idir = 0; idir < 3 ; idir++) ff += 2*n*n*kvec[idir]*kvec[idir]*cos(2*n*dot(kvec, rr));
	return ff;
}

auto gradient_sqwave(inq::math::vector3<double> rr, int n){
	auto kvec = 2.0*M_PI*inq::math::vector3<double>(1.0/9.0, 1.0/12.0, 1.0/10.0);
	inq::math::vector3<double> ff;
	for(int idir = 0; idir < 3 ; idir++) ff[idir] = 2*n*kvec[idir]*cos(n*dot(kvec, rr))*sin(n*dot(kvec, rr));
	return ff;
}

auto gaussiansigma(inq::math::vector3<double> rr, double sigma){
	return  pow(2*M_PI, -1.5)*pow(sigma, -3.0)*exp(-0.5*norm(rr)*pow(sigma, -2.0));
}

auto gaussian(inq::math::vector3<double> rr){  // sigma = 1/sqrt(2)
	return  pow(M_PI, -1.5)*exp(-norm(rr));
}

auto dgaussian(inq::math::vector3<double> rr){
	inq::math::vector3<double> ff;
	for(int idir = 0; idir < 3 ; idir++) ff[idir] = -2.0*rr[idir]*gaussian(rr);
	return ff;
}

auto laplacian_gaussian(inq::math::vector3<double> rr){
	return 4.0*dot(rr, rr)*gaussian(rr) - 6.0*gaussian(rr);
}

TEST_CASE("function hamiltonian::xc_functional", "[hamiltonian::xc_functional]") {

	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
	using Catch::Approx;

	using namespace operations;
	using namespace inq::utils;
	using math::vector3;

	double lx = 9;
	double ly = 12;
	double lz = 10;

	parallel::cartesian_communicator<2> cart_comm(boost::mpi3::environment::get_world_instance(), {});

	SECTION("LDA"){
		systems::box box = systems::box::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b).cutoff_energy(20.0_Ha);
		basis::real_space rs(box, cart_comm);

		basis::field_set<basis::real_space, double> gaussian_field(rs, 1);
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					gaussian_field.hypercubic()[ix][iy][iz][0] = gaussian(vec);
				}
			}
		}

		inq::hamiltonian::xc_functional ldafunctional(XC_LDA_X, 1);
		basis::field_set<basis::real_space, double> gaussianVxc(rs, 1);
		double gaussianExc;

		CHECK(ldafunctional.exx_coefficient() == 0.0);
		
		ldafunctional(gaussian_field, gaussianExc, gaussianVxc);
		CHECK(gaussianExc == -0.270646_a);

		double int_xc_energy = 0.0;
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					math::array<double, 1> local_density{gaussian(vec)};
					math::array<double, 1> local_exc{1};
					math::array<double, 1> local_vxc{1};
					xc_lda_exc_vxc(ldafunctional.libxc_func_ptr(), 1, raw_pointer_cast(local_density.data_elements()), raw_pointer_cast(local_exc.data_elements()), raw_pointer_cast(local_vxc.data_elements()));
					gpu::sync();
					
					CHECK(Approx(local_vxc[0]) == gaussianVxc.hypercubic()[ix][iy][iz][0]);

					int_xc_energy += local_exc[0]*local_density[0]*rs.volume_element();
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&int_xc_energy, 1, std::plus<>{});
		
		CHECK(Approx(gaussianExc) == int_xc_energy);

		if(rs.part().contains(1)) CHECK(gaussianVxc.matrix()[rs.part().global_to_local(parallel::global_index(1))][0] == -0.5111609291_a);
		if(rs.part().contains(8233)) CHECK(gaussianVxc.matrix()[rs.part().global_to_local(parallel::global_index(8233))][0] == -0.00406881_a);
		if(rs.part().contains(233)) CHECK(fabs(gaussianVxc.matrix()[rs.part().global_to_local(parallel::global_index(233))][0]) < 1e-10);
		if(rs.part().contains(rs.size() - 1)) CHECK(gaussianVxc.matrix()[rs.part().global_to_local(parallel::global_index(rs.size() - 1))][0] == -0.4326883849_a);

	}

	SECTION("LSDA"){
		
		systems::box box = systems::box::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b).cutoff_energy(20.0_Ha);
		basis::real_space rs(box, cart_comm);

		basis::field_set<basis::real_space, double> density_unp(rs, 2);
		basis::field_set<basis::real_space, double> density_pol(rs, 2);		
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					density_unp.hypercubic()[ix][iy][iz][0] = 0.5*gaussian(vec);
					density_unp.hypercubic()[ix][iy][iz][1] = 0.5*gaussian(vec);
					density_pol.hypercubic()[ix][iy][iz][0] = 0.3*gaussian(vec);
					density_pol.hypercubic()[ix][iy][iz][1] = 0.7*gaussian(vec);					
				}
			}
		}

		inq::hamiltonian::xc_functional ldafunctional(XC_LDA_X, 2);

		basis::field_set<basis::real_space, double> vxc_unp(rs, 2);
		basis::field_set<basis::real_space, double> vxc_pol(rs, 2);		

		double exc_unp, exc_pol;

		CHECK(ldafunctional.exx_coefficient() == 0.0);
		
		ldafunctional(density_unp, exc_unp, vxc_unp);
		ldafunctional(density_pol, exc_pol, vxc_pol);		

		CHECK(exc_unp == -0.270646_a);
		CHECK(exc_pol == -0.2804198447_a);		
		
		double int_exc_unp = 0.0;
		double int_exc_pol = 0.0;		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					auto vec = rs.point_op().rvector_cartesian(ix, iy, iz);
					math::array<double, 1> local_density{0.5*gaussian(vec), 0.5*gaussian(vec)};
					math::array<double, 1> local_exc{NAN};
					math::array<double, 1> local_vxc{NAN, NAN};
					
					xc_lda_exc_vxc(ldafunctional.libxc_func_ptr(), 1, raw_pointer_cast(local_density.data_elements()), raw_pointer_cast(local_exc.data_elements()), raw_pointer_cast(local_vxc.data_elements()));
					gpu::sync();

					CHECK(vxc_unp.hypercubic()[ix][iy][iz][0] == Approx(vxc_unp.hypercubic()[ix][iy][iz][1]));					
					CHECK(Approx(local_vxc[0]) == vxc_unp.hypercubic()[ix][iy][iz][0]);
					CHECK(Approx(local_vxc[1]) == vxc_unp.hypercubic()[ix][iy][iz][1]);

					int_exc_unp += local_exc[0]*(local_density[0] + local_density[1])*rs.volume_element();
					
					local_density[0] = 0.7*gaussian(vec);
					local_density[1] = 0.3*gaussian(vec);					
					xc_lda_exc_vxc(ldafunctional.libxc_func_ptr(), 1, raw_pointer_cast(local_density.data_elements()), raw_pointer_cast(local_exc.data_elements()), raw_pointer_cast(local_vxc.data_elements()));
					gpu::sync();

					CHECK(fabs(local_vxc[1] - vxc_pol.hypercubic()[ix][iy][iz][0]) < 1e-14);
					CHECK(fabs(local_vxc[0] - vxc_pol.hypercubic()[ix][iy][iz][1]) < 1e-14);

					int_exc_pol += local_exc[0]*(local_density[0] + local_density[1])*rs.volume_element();
				}
			}
		}

		cart_comm.all_reduce_in_place_n(&int_exc_unp, 1, std::plus<>{});
		cart_comm.all_reduce_in_place_n(&int_exc_pol, 1, std::plus<>{});

		CHECK(Approx(exc_unp) == int_exc_unp);
		CHECK(Approx(exc_pol) == int_exc_pol);

		if(rs.part().contains(1)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(1))][0] == -0.5111609291_a);
		if(rs.part().contains(1)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(1))][1] == -0.5111609291_a);		
		if(rs.part().contains(8233)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(8233))][0] == -0.00406881_a);
		if(rs.part().contains(8233)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(8233))][1] == -0.00406881_a);		
		if(rs.part().contains(233)) CHECK(fabs(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(233))][0]) < 1e-10);
		if(rs.part().contains(233)) CHECK(fabs(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(233))][1]) < 1e-10);		
		if(rs.part().contains(rs.size() - 1)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(rs.size() - 1))][0] == -0.4326883849_a);
		if(rs.part().contains(rs.size() - 1)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(rs.size() - 1))][1] == -0.4326883849_a);		

		if(rs.part().contains(1)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(1))][0] == -0.4311298248_a);
		if(rs.part().contains(1)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(1))][1] == -0.571830079_a);		
		if(rs.part().contains(8233)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(8233))][0] == -0.0034317692_a);
		if(rs.part().contains(8233)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(8233))][1] == -0.0045517353_a);		
		if(rs.part().contains(233)) CHECK(fabs(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(233))][0]) < 1e-10);
		if(rs.part().contains(233)) CHECK(fabs(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(233))][1]) < 1e-10);		
		if(rs.part().contains(rs.size() - 1)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(rs.size() - 1))][0] == -0.3649435178_a);
		if(rs.part().contains(rs.size() - 1)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(rs.size() - 1))][1] == -0.4840437116_a);
	}
		
	SECTION("GGA"){
		systems::box box = systems::box::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b).cutoff_energy(90.0_Ha);
		basis::real_space rs(box, cart_comm);
		
		basis::field_set<basis::real_space, double> field(rs, 1);
		gpu::run(rs.local_sizes()[2], rs.local_sizes()[1], rs.local_sizes()[0],
						 [pop = rs.point_op(), fie = begin(field.hypercubic())] GPU_LAMBDA (auto iz, auto iy, auto ix){
							 auto vec = pop.rvector_cartesian(ix, iy, iz);
							 fie[ix][iy][iz][0] = sqwave(vec, 3);
						 });
	
		inq::hamiltonian::xc_functional ggafunctional(XC_GGA_X_PBE, 1);
		basis::field_set<basis::real_space, double> Vxc(rs, 1);

		CHECK(ggafunctional.exx_coefficient() == 0.0);
		
		double Exc = 0.0;
		ggafunctional(field, Exc, Vxc);

		CHECK(Exc == -393.4604748792_a);
		
		if(rs.part().contains(1)) CHECK(Vxc.matrix()[rs.part().global_to_local(parallel::global_index(1))][0] == -0.5607887985_a);
		if(rs.part().contains(33)) CHECK(Vxc.matrix()[rs.part().global_to_local(parallel::global_index(33))][0] == -1.1329131862_a);
		if(rs.part().contains(rs.size() - 1)) CHECK(Vxc.matrix()[rs.part().global_to_local(parallel::global_index(rs.size() - 1))][0] == -1.1461742979_a);
	}

	SECTION("Spin GGA"){
		
		systems::box box = systems::box::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b).cutoff_energy(90.0_Ha);
		basis::real_space rs(box, cart_comm);
		
		basis::field_set<basis::real_space, double> density_unp(rs, 2);
		basis::field_set<basis::real_space, double> density_pol(rs, 2);		
		
		gpu::run(rs.local_sizes()[2], rs.local_sizes()[1], rs.local_sizes()[0],
						 [pop = rs.point_op(), den_unp = begin(density_unp.hypercubic()), den_pol = begin(density_pol.hypercubic())] GPU_LAMBDA (auto iz, auto iy, auto ix){
							 auto vec = pop.rvector_cartesian(ix, iy, iz);
							 den_unp[ix][iy][iz][0] = 0.5*sqwave(vec, 3);
							 den_unp[ix][iy][iz][1] = 0.5*sqwave(vec, 3);
							 den_pol[ix][iy][iz][0] = 0.3*sqwave(vec, 3);
							 den_pol[ix][iy][iz][1] = 0.7*sqwave(vec, 3);							 
						 });
	
		inq::hamiltonian::xc_functional ggafunctional(XC_GGA_X_PBE, 2);
		
		basis::field_set<basis::real_space, double> vxc_unp(rs, 2);
		basis::field_set<basis::real_space, double> vxc_pol(rs, 2);

		CHECK(ggafunctional.exx_coefficient() == 0.0);
		
		double exc_unp, exc_pol;
		
		ggafunctional(density_unp, exc_unp, vxc_unp);
		ggafunctional(density_pol, exc_pol, vxc_pol);		

		CHECK(exc_unp == -393.4604748792_a);
		CHECK(exc_pol == -406.0224635643_a);
		
		if(rs.part().contains(1)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(1))][0] == -0.5607887985_a);
		if(rs.part().contains(1)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(1))][1] == -0.5607887985_a);		
		if(rs.part().contains(33)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(33))][0] == -1.1329131862_a);
		if(rs.part().contains(33)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(33))][1] == -1.1329131862_a);		
		if(rs.part().contains(rs.size() - 1)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(rs.size() - 1))][0] == -1.1461742979_a);
		if(rs.part().contains(rs.size() - 1)) CHECK(vxc_unp.matrix()[rs.part().global_to_local(parallel::global_index(rs.size() - 1))][1] == -1.1461742979_a);

		if(rs.part().contains(1)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(1))][0] == -0.4941651943_a);
		if(rs.part().contains(1)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(1))][1] == -0.6129310854_a);		
		if(rs.part().contains(33)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(33))][0] == -0.9675738833_a);
		if(rs.part().contains(33)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(33))][1] == -1.2504476169_a);		
		if(rs.part().contains(rs.size() - 1)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(rs.size() - 1))][0] == -1.0006730841_a);
		if(rs.part().contains(rs.size() - 1)) CHECK(vxc_pol.matrix()[rs.part().global_to_local(parallel::global_index(rs.size() - 1))][1] == -1.252993172_a);		
	}
	
	SECTION("Uniform"){
		systems::box box = systems::box::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b).cutoff_energy(20.0_Ha);
		basis::real_space rs(box, cart_comm);
		
		basis::field_set<basis::real_space, double> gaussian_field(rs, 1);
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					gaussian_field.hypercubic()[ix][iy][iz][0] = 0.393;
				}
			}
		}
	
		inq::hamiltonian::xc_functional ggafunctional(XC_GGA_X_PBE, 1);
		inq::hamiltonian::xc_functional ldafunctional(XC_LDA_X, 1);
		basis::field_set<basis::real_space, double> gaussianVxcLDA(rs, 1), gaussianVxcGGA(rs, 1);
		double gaussianExcLDA, gaussianExcGGA;

		ggafunctional(gaussian_field, gaussianExcLDA, gaussianVxcLDA);
		ldafunctional(gaussian_field, gaussianExcGGA, gaussianVxcGGA);
		CHECK(gaussianExcLDA == gaussianExcGGA);
		
		for(int ix = 0; ix < rs.local_sizes()[0]; ix++){
			for(int iy = 0; iy < rs.local_sizes()[1]; iy++){
				for(int iz = 0; iz < rs.local_sizes()[2]; iz++){
					CHECK(Approx(gaussianVxcLDA.hypercubic()[ix][iy][iz][0]) == gaussianVxcGGA.hypercubic()[ix][iy][iz][0]);
				}
			}
		}

	}

	inq::hamiltonian::xc_functional b3lyp(XC_HYB_GGA_XC_B3LYP, 1);
	inq::hamiltonian::xc_functional pbeh(XC_HYB_GGA_XC_PBEH, 1);
	
	SECTION("HYBRIDS"){
		CHECK(b3lyp.exx_coefficient() == 0.2_a);
		CHECK(pbeh.exx_coefficient() == 0.25_a);
	}

	SECTION("COPY AND ASSIGNMENT"){
		auto copy = b3lyp;
		CHECK(copy.exx_coefficient() == 0.2_a);

		copy = pbeh;
		CHECK(copy.exx_coefficient() == 0.25_a);

		auto copy2 = std::move(copy);
		CHECK(copy2.exx_coefficient() == 0.25_a);

		copy2 = std::move(b3lyp);
		CHECK(copy2.exx_coefficient() == 0.2_a);
		
	}
	
}

#endif
#endif
