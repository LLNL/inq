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
#include <observables/magnetization.hpp>
#include <operations/add.hpp>
#include <operations/integral.hpp>
#include <options/theory.hpp>
#include <hamiltonian/xc_functional.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <perturbations/none.hpp>
#include <utils/profiling.hpp>
#include <systems/electrons.hpp>
#include <cmath>

namespace inq {
namespace hamiltonian {

class xc_term {

	std::vector<hamiltonian::xc_functional> functionals_;	

public:
	
	xc_term(options::theory interaction, int const spin_components){
		functionals_.emplace_back(int(interaction.exchange()), spin_components);
		functionals_.emplace_back(int(interaction.correlation()), spin_components);
	}

  ////////////////////////////////////////////////////////////////////////////////////////////

	auto any_requires_gradient() const {
		for(auto & func : functionals_) if(func.requires_gradient()) return true;
		return false;
	}

	////////////////////////////////////////////////////////////////////////////////////////////

	auto any_true_functional() const {
		for(auto & func : functionals_) if(func.true_functional()) return true;
		return false;
	}
	
  ////////////////////////////////////////////////////////////////////////////////////////////

	template <typename SpinDensityType, typename CoreDensityType>
	SpinDensityType process_density(SpinDensityType const & spin_density, CoreDensityType const & core_density) const{

		SpinDensityType full_density(spin_density.basis(), std::min(2, spin_density.set_size()));

		if(spin_density.set_size() == 4) {
			gpu::run(spin_density.basis().local_size(),
							 [spi = begin(spin_density.matrix()), ful = begin(full_density.matrix()), cor = begin(core_density.linear())] GPU_LAMBDA (auto ip){
								 auto dtot = spi[ip][0] + spi[ip][1];
								 auto mag = observables::local_magnetization(spi[ip], 4);
								 auto dpol = mag.length();
								 ful[ip][0] = 0.5*(dtot + dpol);
								 ful[ip][1] = 0.5*(dtot - dpol);
								 for(int ispin = 0; ispin < 2; ispin++){
									 if(ful[ip][ispin] < 0.0) ful[ip][ispin] = 0.0;
									 ful[ip][ispin] += cor[ip]/2;
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

	template <typename SpinDensityType, typename VXC, typename VKS>
	void process_potential(SpinDensityType const & spin_density, VXC const & vxc, VKS & vks) const {
		
		if (spin_density.set_size() == 4) {
			gpu::run(vxc.basis().local_size(),
						 [spi = begin(spin_density.matrix()), vx = begin(vxc.matrix()), vk = begin(vks.matrix())] GPU_LAMBDA (auto ip){
							auto bxc = 0.5*(vx[ip][0] - vx[ip][1]);
							auto vxc = 0.5*(vx[ip][0] + vx[ip][1]);
							auto mag = observables::local_magnetization(spi[ip], 4);
							auto dpol = mag.length();
							if (fabs(dpol) > 1.e-7) {
								auto e_mag = mag/dpol;
								vk[ip][0] += vxc + bxc*e_mag[2];
								vk[ip][1] += vxc - bxc*e_mag[2];
								vk[ip][2] += bxc*e_mag[0];
								vk[ip][3] += bxc*e_mag[1];
							}
							else {
								vk[ip][0] += vxc;
								vk[ip][1] += vxc;
							}
						 });
		}
		else {
			assert(spin_density.set_size() == 1 or spin_density.set_size() == 2);
			gpu::run(vxc.local_set_size(), vxc.basis().local_size(),
						 [vx = begin(vxc.matrix()), vk = begin(vks.matrix())] GPU_LAMBDA (auto is, auto ip){
							vk[ip][is] += vx[ip][is];
						 });
		}

	}

  ///////////////////////////////////////////////////////////////////////////////////////////

	template <typename SpinDensityType, typename VXC>
	auto compute_nvxc(SpinDensityType const & spin_density, VXC const & vfunc) const {

		auto nvxc_ = 0.;
		if (spin_density.set_size() == 4) {
			basis::field_set<basis::real_space, double> vxc(spin_density.skeleton());
			vxc.fill(0.0);
			gpu::run(vfunc.basis().local_size(),
						 [spi = begin(spin_density.matrix()), vxi = begin(vfunc.matrix()), vxf = begin(vxc.matrix())] GPU_LAMBDA (auto ip){
							auto b = 0.5*(vxi[ip][0] - vxi[ip][1]);
							auto v = 0.5*(vxi[ip][0] + vxi[ip][1]);
							auto mag = observables::local_magnetization(spi[ip], 4);
							auto dpol = mag.length();
							if (fabs(dpol) > 1.e-7) {
								auto e_mag = mag/dpol;
								// Vxc = [vxc+bxc^z, bxc^-; bxc^+, vxc-bxc^z]
								vxf[ip][0] = v + b*e_mag[2];
								vxf[ip][1] = v - b*e_mag[2];
								vxf[ip][2] = 2.0*b*e_mag[0];
								vxf[ip][3] = 2.0*b*e_mag[1];
							}
							else {
								vxf[ip][0] = v;
								vxf[ip][1] = v;
							}
						 });
			basis::field<basis::real_space, double> rfield(spin_density.basis());
			rfield.fill(0.0);
			gpu::run(spin_density.local_set_size(), spin_density.basis().local_size(),
						[spi = begin(spin_density.matrix()), vx = begin(vxc.matrix()), rf = begin(rfield.linear())] GPU_LAMBDA (auto is, auto ip){
							rf[ip] += vx[ip][is] * spi[ip][is];
						});
			nvxc_ += operations::integral(rfield);
		}
		else {
			nvxc_ = operations::integral_product_sum(spin_density, vfunc);
		}
		return nvxc_;
	}

  ////////////////////////////////////////////////////////////////////////////////////////////
	
  template <typename SpinDensityType, typename CoreDensityType, typename VKSType>
  void operator()(SpinDensityType const & spin_density, CoreDensityType const & core_density, VKSType & vks, double & exc, double & nvxc) const {
    
    exc = 0.0;
		nvxc = 0.0;
		if(not any_true_functional()) return;
		
		auto full_density = process_density(spin_density, core_density);
		
		double efunc = 0.0;
		
		basis::field_set<basis::real_space, double> vfunc(full_density.skeleton());
		auto density_gradient = std::optional<decltype(operations::gradient(full_density))>{};
		if(any_requires_gradient()) density_gradient.emplace(operations::gradient(full_density));

		for(auto & func : functionals_){
			if(not func.true_functional()) continue;

			evaluate_functional(func, full_density, density_gradient, efunc, vfunc);
			exc += efunc;
			process_potential(spin_density, vfunc, vks);

			nvxc += compute_nvxc(spin_density, vfunc);
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////

	template <typename SpinDensityType, typename CoreDensityType, typename VxcType>
	void compute_vxc(SpinDensityType const & spin_density, CoreDensityType const & core_density, VxcType & vxc) {

		if(not any_true_functional()) return;
		auto full_density = process_density(spin_density, core_density);
		double efunc = 0.0;
		basis::field_set<basis::real_space, double> vfunc(full_density.skeleton());
		vfunc.fill(0.0);
		vxc.fill(0.0);
		auto density_gradient = std::optional<decltype(operations::gradient(full_density))>{};
		if(any_requires_gradient()) density_gradient.emplace(operations::gradient(full_density));

		for(auto & func : functionals_) {
			if(not func.true_functional()) continue;
			evaluate_functional(func, full_density, density_gradient, efunc, vfunc);
			if (spin_density.set_size() == 4) {
				gpu::run(vfunc.basis().local_size(),
				[spi = begin(spin_density.matrix()), vxi = begin(vfunc.matrix()), vxf = begin(vxc.matrix())] GPU_LAMBDA (auto ip){
					auto b = 0.5*(vxi[ip][0] - vxi[ip][1]);
					auto v = 0.5*(vxi[ip][0] + vxi[ip][1]);
					auto mag = observables::local_magnetization(spi[ip], 4);
					auto dpol = mag.length();
					if (fabs(dpol) > 1.e-7) {
						auto e_mag = mag / dpol;
						vxf[ip][0] += v + b*e_mag[2];
						vxf[ip][1] += v - b*e_mag[2];
						vxf[ip][2] += b*e_mag[0];
						vxf[ip][3] += b*e_mag[1];
					}
					else {
						vxf[ip][0] += v;
						vxf[ip][1] += v;
					}
				});
			}
			else {
				assert(spin_density.set_size() == 1 or spin_density.set_size() == 2);
				gpu::run(vfunc.local_set_size(), vfunc.basis().local_size(),
					[vxi = begin(vfunc.matrix()), vxf = begin(vxc.matrix())] GPU_LAMBDA (auto is, auto ip){
						vxf[ip][is] += vxi[ip][is];
					});
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////

	template <typename DensityType, typename DensityGradientType>
	static void evaluate_functional(hamiltonian::xc_functional const & functional, DensityType const & density, DensityGradientType const & density_gradient,
																	double & efunctional, basis::field_set<basis::real_space, double> & vfunctional){
		CALI_CXX_MARK_FUNCTION;

		auto edens = basis::field<basis::real_space, double>(density.basis());

		assert(functional.nspin() == density.set_size());
		
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
		return functionals_[0];
	}
	
  ////////////////////////////////////////////////////////////////////////////////////////////
	
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

	if(comm.size() > 4) return; //FIXME: check the problem for size 5. It returns a multi error
	
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
			CHECK(vfunc_unp.matrix()[index][0] == 0.0005467193_a);
			CHECK(vfunc_pol.matrix()[index][0] == 0.0005956583_a);
			CHECK(vfunc_pol.matrix()[index][1] == 0.0005978958_a);
		}

		if(bas.part().contains(4444)) {
			auto index = bas.part().global_to_local(parallel::global_index(4444));
			CHECK(vfunc_unp.matrix()[index][0] == -0.0798456253_a);
			CHECK(vfunc_pol.matrix()[index][0] == -0.0667968142_a);
			CHECK(vfunc_pol.matrix()[index][1] == -0.0830118308_a);
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
			CHECK(vfunc_unp.matrix()[index][0] == -0.6495909727_a);
			CHECK(vfunc_pol.matrix()[index][0] == -0.6398010386_a);
			CHECK(vfunc_pol.matrix()[index][1] == -0.6142058762_a);
		}

		if(bas.part().contains(4444)) {
			auto index = bas.part().global_to_local(parallel::global_index(4444));
			CHECK(vfunc_unp.matrix()[index][0] == -0.2879332051_a);
			CHECK(vfunc_pol.matrix()[index][0] == -0.3195127242_a);
			CHECK(vfunc_pol.matrix()[index][1] == -0.2368583776_a);
		}

	}

	SECTION("xc_term object") {

		auto hf = hamiltonian::xc_term(options::theory{}.hartree_fock(), 1);
		CHECK(hf.any_requires_gradient() == false);
		CHECK(hf.any_true_functional() == false);
		
		auto lda = hamiltonian::xc_term(options::theory{}.lda(), 1);
		CHECK(lda.any_requires_gradient() == false);
		CHECK(lda.any_true_functional() == true);

		auto pbe = hamiltonian::xc_term(options::theory{}.pbe(), 1);
		CHECK(pbe.any_requires_gradient() == true);
		CHECK(pbe.any_true_functional() == true);

	}
	
}
#endif
