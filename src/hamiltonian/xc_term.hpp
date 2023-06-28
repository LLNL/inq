/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__XC_TERM
#define INQ__HAMILTONIAN__XC_TERM

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <basis/field.hpp>
#include <solvers/poisson.hpp>
#include <observables/density.hpp>
#include <operations/add.hpp>
#include <operations/integral.hpp>
#include <options/theory.hpp>
#include <hamiltonian/xc_functional.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <perturbations/none.hpp>
#include <utils/profiling.hpp>

namespace inq {
namespace hamiltonian {

class xc_term {
	
public:
	
	xc_term(options::theory interaction, int const spin_components):
		exchange_(int(interaction.exchange()), spin_components),
		correlation_(int(interaction.correlation()), spin_components)
	{
	}

  ////////////////////////////////////////////////////////////////////////////////////////////

	template <typename SpinDensityType, typename CoreDensityType>
	SpinDensityType process_density(SpinDensityType const & spin_density, CoreDensityType const & core_density) const{

		SpinDensityType full_density(spin_density.basis(), std::min(2, spin_density.set_size()));

		if(spin_density.set_size() == 4) {
			//for spinors convert the density to 2 components
			gpu::run(spin_density.basis().local_size(),
							 [spi = begin(spin_density.matrix()), ful = begin(full_density.matrix()), cor = begin(core_density.linear())] GPU_LAMBDA (auto ip){
								 auto dtot = spi[ip][0] + spi[ip][1];
								 auto dd = spi[ip][0] - spi[ip][1];
								 auto dpol = sqrt(dd*dd + 4.0*(spi[ip][2]*spi[ip][2] + spi[ip][3]*spi[ip][3]));
								 ful[ip][0] = 0.5*(dtot + dpol);
								 ful[ip][1] = 0.5*(dtot - dpol);
								 for(int ispin = 0; ispin < 2; ispin++){
									 if(ful[ip][ispin] < 0.0) ful[ip][ispin] = 0.0;
									 ful[ip][ispin] += 0.5*cor[ip];
								 }
							 });
		} else {
			assert(spin_density.set_size() == 1 or spin_density.set_size() == 2);
			
			gpu::run(spin_density.basis().local_size(),
							 [spi = begin(spin_density.matrix()), ful = begin(full_density.matrix()), cor = begin(core_density.linear()), nspin = spin_density.set_size()] GPU_LAMBDA (auto ip){
								 for(int ispin = 0; ispin < nspin; ispin++){
									 ful[ip][ispin] = spi[ip][ispin];
									 if(ful[ip][ispin] < 0.0) ful[ip][ispin] = 0.0;
									 ful[ip][ispin] += cor[ip]/nspin;
								 }
							 });
		}
		
		return full_density;
	}
	
  ////////////////////////////////////////////////////////////////////////////////////////////
	
  template <typename SpinDensityType, typename CoreDensityType, typename VKSType>
  void operator()(SpinDensityType const & spin_density, CoreDensityType const & core_density, VKSType & vks, double & exc, double & nvxc) const {
    
    exc = 0.0;
		nvxc = 0.0;
		if(not exchange_.true_functional() and not correlation_.true_functional()) return;
		
		auto full_density = process_density(spin_density, core_density);
		
		double efunc = 0.0;
		basis::field_set<basis::real_space, double> vfunc(spin_density.skeleton());

		auto density_gradient = std::optional<decltype(operations::gradient(full_density))>{};
		if(exchange_.requires_gradient() or correlation_.requires_gradient()){
			density_gradient.emplace(operations::gradient(full_density));
		}
			
		if(exchange_.true_functional()){
			evaluate_functional(exchange_, full_density, density_gradient, efunc, vfunc);
			exc += efunc;
			operations::increment(vks, vfunc);
			nvxc += operations::integral_product_sum(spin_density, vfunc); //the core correction does not go here
		}
		
		if(correlation_.true_functional()){
			evaluate_functional(correlation_, full_density, density_gradient, efunc, vfunc);
			exc += efunc;
			operations::increment(vks, vfunc);
			nvxc += operations::integral_product_sum(spin_density, vfunc); //the core correction does not go here
		}
  }

	////////////////////////////////////////////////////////////////////////////////////////////

	template <typename DensityType, typename DensityGradientType>
	static void evaluate_functional(hamiltonian::xc_functional const & functional, DensityType const & density, DensityGradientType const & density_gradient,
													 double & efunctional, basis::field_set<basis::real_space, double> & vfunctional){
		CALI_CXX_MARK_FUNCTION;

		auto edens = basis::field<basis::real_space, double>(density.basis());
		
		if(functional.family() == XC_FAMILY_LDA){
			
			xc_lda_exc_vxc(functional.libxc_func_ptr(), density.basis().local_size(), raw_pointer_cast(density.matrix().data_elements()),
										 raw_pointer_cast(edens.linear().data_elements()), raw_pointer_cast(vfunctional.matrix().data_elements()));
			gpu::sync();
			
		} else if(functional.family() == XC_FAMILY_GGA){

			auto nsig = (density.set_size() > 1) ? 3:1;
			
			basis::field_set<basis::real_space, double> sig(density.basis(), nsig);
			basis::field_set<basis::real_space, double> vsig(sig.skeleton());

			gpu::run(density.basis().local_size(),
							 [gr = begin(density_gradient->matrix()), si = begin(sig.matrix()), metric = density.basis().cell().metric(), nsig] GPU_LAMBDA (auto ip){
								 si[ip][0] = metric.norm(gr[ip][0]);
								 if(nsig > 1) si[ip][1] = metric.dot(gr[ip][0], gr[ip][1]);
								 if(nsig > 1) si[ip][2] = metric.norm(gr[ip][1]);
							 });

			xc_gga_exc_vxc(functional.libxc_func_ptr(), density.basis().local_size(), raw_pointer_cast(density.matrix().data_elements()), raw_pointer_cast(sig.matrix().data_elements()),
										 raw_pointer_cast(edens.linear().data_elements()), raw_pointer_cast(vfunctional.matrix().data_elements()), raw_pointer_cast(vsig.matrix().data_elements()));
			gpu::sync();

			basis::field_set<basis::real_space, vector3<double, covariant>> term(vfunctional.skeleton());

			gpu::run(density.basis().local_size(),
							 [vs = begin(vsig.matrix()), gr = begin(density_gradient->matrix()), te = begin(term.matrix()), nsig] GPU_LAMBDA (auto ip){
								 if(nsig == 1) te[ip][0] = -2.0*vs[ip][0]*gr[ip][0];
								 if(nsig == 3) te[ip][0] = -2.0*vs[ip][0]*gr[ip][0] - vs[ip][1]*gr[ip][1];
								 if(nsig == 3) te[ip][1] = -2.0*vs[ip][2]*gr[ip][1] - vs[ip][1]*gr[ip][0];
							 });

			auto div_term = operations::divergence(term);

			gpu::run(density.local_set_size(), density.basis().local_size(),
							 [di = begin(div_term.matrix()), vf = begin(vfunctional.matrix())] GPU_LAMBDA (auto ispin, auto ip){
								 vf[ip][ispin] += di[ip][ispin];
							 });

		} else {
			std::runtime_error("inq error: unsupported exchange correlation functional type");
		}
		
		efunctional = operations::integral_product(edens, observables::density::total(density));
		
	}
	
  ////////////////////////////////////////////////////////////////////////////////////////////
	
	auto & exchange() const {
		return exchange_;
	}
	
  ////////////////////////////////////////////////////////////////////////////////////////////
	
private:
	hamiltonian::xc_functional exchange_;
	hamiltonian::xc_functional correlation_;
	
};
}
}
#endif

#ifdef INQ_HAMILTONIAN_XC_TERM_UNIT_TEST
#undef INQ_HAMILTONIAN_XC_TERM_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <basis/real_space.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG){

	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
	using namespace operations;
	using Catch::Approx;
	
	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

	auto lx = 10.3;
	auto ly = 13.8;
	auto lz =  4.5;
	
	basis::real_space bas(systems::cell::orthorhombic(lx*1.0_b, ly*1.0_b, lz*1.0_b), /*spacing =*/ 0.40557787, comm);

	basis::field_set<basis::real_space, double> density_unp(bas, 1);	
	basis::field_set<basis::real_space, double> density_pol(bas, 2);
	
	//Define k-vector for test function
	auto kvec = 2.0*M_PI*vector3<double>(1.0/lx, 1.0/ly, 1.0/lz);
	
	auto ff = [] (auto & kk, auto & rr){
		return std::max(0.0, cos(dot(kk, rr)) + 1.0);
	};

	for(int ix = 0; ix < bas.local_sizes()[0]; ix++){
		for(int iy = 0; iy < bas.local_sizes()[1]; iy++){
			for(int iz = 0; iz < bas.local_sizes()[2]; iz++){
				auto vec = bas.point_op().rvector_cartesian(ix, iy, iz);
				density_unp.hypercubic()[ix][iy][iz][0] = ff(kvec, vec);
				auto pol = sin(norm(vec)/100.0);
				density_pol.hypercubic()[ix][iy][iz][0] = (1.0 - pol)*ff(kvec, vec);
				density_pol.hypercubic()[ix][iy][iz][1] = pol*ff(kvec, vec);
			}
		}
	}

	observables::density::normalize(density_unp, 42.0);
	observables::density::normalize(density_pol, 42.0);
	
	CHECK(operations::integral_sum(density_unp) == 42.0_a);
	CHECK(operations::integral_sum(density_pol) == 42.0_a);	
	
	auto grad_unp = std::optional{operations::gradient(density_unp)};
	auto grad_pol = std::optional{operations::gradient(density_pol)};

	if(bas.part().contains(5439)) {
		auto index = bas.part().global_to_local(parallel::global_index(5439));
		CHECK(density_unp.matrix()[index][0] == 0.0024885602_a);
		CHECK(density_pol.matrix()[index][0] == 0.0009452194_a);
		CHECK(density_pol.matrix()[index][1] == 0.0015433408_a);
	}

	basis::field_set<basis::real_space, double> vfunc_unp(bas, 1);	
	basis::field_set<basis::real_space, double> vfunc_pol(bas, 2);
	
	SECTION("LDA_X"){
		
		hamiltonian::xc_functional func_unp(XC_LDA_X, 1);
		hamiltonian::xc_functional func_pol(XC_LDA_X, 2);
		
		double efunc_unp = NAN;
		double efunc_pol = NAN;
		
		hamiltonian::xc_term::evaluate_functional(func_unp, density_unp, grad_unp, efunc_unp, vfunc_unp);
		hamiltonian::xc_term::evaluate_functional(func_pol, density_pol, grad_pol, efunc_pol, vfunc_pol);

		CHECK(efunc_unp == -14.0558385758_a);
		CHECK(efunc_pol == -15.1704508993_a);

		if(bas.part().contains(5439)) {
			auto index = bas.part().global_to_local(parallel::global_index(5439));
			CHECK(vfunc_unp.matrix()[index][0] == -0.1334462916_a);
			CHECK(vfunc_pol.matrix()[index][0] == -0.1217618773_a);
			CHECK(vfunc_pol.matrix()[index][1] == -0.1433797225_a);
		}

		if(bas.part().contains(4444)) {
			auto index = bas.part().global_to_local(parallel::global_index(4444));
			CHECK(vfunc_unp.matrix()[index][0] == -0.3276348215_a);
			CHECK(vfunc_pol.matrix()[index][0] == -0.3784052378_a);
			CHECK(vfunc_pol.matrix()[index][1] == -0.2527984139_a);
		}

	}

	SECTION("PBE_C"){
		
		hamiltonian::xc_functional func_unp(XC_GGA_C_PBE, 1);
		hamiltonian::xc_functional func_pol(XC_GGA_C_PBE, 2);
		
		double efunc_unp = NAN;
		double efunc_pol = NAN;
		
		hamiltonian::xc_term::evaluate_functional(func_unp, density_unp, grad_unp, efunc_unp, vfunc_unp);
		hamiltonian::xc_term::evaluate_functional(func_pol, density_pol, grad_pol, efunc_pol, vfunc_pol);

		CHECK(efunc_unp == -1.8220292936_a);
		CHECK(efunc_pol == -1.5664843681_a);

		if(bas.part().contains(5439)) {
			auto index = bas.part().global_to_local(parallel::global_index(5439));
			CHECK(vfunc_unp.matrix()[index][0] == 0.0005665189_a);
			CHECK(vfunc_pol.matrix()[index][0] == 0.0007334274_a);
			CHECK(vfunc_pol.matrix()[index][1] == 0.0007356649_a);
		}

		if(bas.part().contains(4444)) {
			auto index = bas.part().global_to_local(parallel::global_index(4444));
			CHECK(vfunc_unp.matrix()[index][0] == -0.0770230402_a);
			CHECK(vfunc_pol.matrix()[index][0] == -0.0641968735_a);
			CHECK(vfunc_pol.matrix()[index][1] == -0.0804118901_a);
		}

	}

	SECTION("B3LYP"){
		
		hamiltonian::xc_functional func_unp(XC_HYB_GGA_XC_B3LYP, 1);
		hamiltonian::xc_functional func_pol(XC_HYB_GGA_XC_B3LYP, 2);
		
		double efunc_unp = NAN;
		double efunc_pol = NAN;
		
		hamiltonian::xc_term::evaluate_functional(func_unp, density_unp, grad_unp, efunc_unp, vfunc_unp);
		hamiltonian::xc_term::evaluate_functional(func_pol, density_pol, grad_pol, efunc_pol, vfunc_pol);

		CHECK(efunc_unp == -13.2435562623_a);
		CHECK(efunc_pol == -13.8397387159_a);

		if(bas.part().contains(5439)) {
			auto index = bas.part().global_to_local(parallel::global_index(5439));
			CHECK(vfunc_unp.matrix()[index][0] == -0.5955644446_a);
			CHECK(vfunc_pol.matrix()[index][0] == -0.5721375366_a);
			CHECK(vfunc_pol.matrix()[index][1] == -0.5667836757_a);
		}

		if(bas.part().contains(4444)) {
			auto index = bas.part().global_to_local(parallel::global_index(4444));
			CHECK(vfunc_unp.matrix()[index][0] == -0.2720801967_a);
			CHECK(vfunc_pol.matrix()[index][0] == -0.3121236386_a);
			CHECK(vfunc_pol.matrix()[index][1] == -0.1997021271_a);
		}

	}
}
#endif
