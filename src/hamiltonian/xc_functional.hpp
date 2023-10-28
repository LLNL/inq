/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__XC_FUNCTIONAL
#define INQ__HAMILTONIAN__XC_FUNCTIONAL

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa, Alexey Kartsev
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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

			if(true_functional()) assert(nspin_ == 1 or nspin_ == 2);
			
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
		
		auto & family() const {
			return func_.info->family;
		}

		auto & libxc_func() const {
			return func_;
		}
		
		auto libxc_func_ptr() const {
			return &func_;
		}

		auto requires_gradient() const {
			if(not true_functional()) return false;
			return family() != XC_FAMILY_LDA;
		}
		
		auto exx_coefficient() const {
			if(xc_hyb_type(libxc_func_ptr()) == XC_HYB_HYBRID) return xc_hyb_exx_coef(libxc_func_ptr());
			return 0.0;
		}

	private:
		
		int id_;
		int nspin_;
		xc_func_type func_;
	};

}
}
#endif

#ifdef INQ_HAMILTONIAN_XC_FUNCTIONAL_UNIT_TEST
#undef INQ_HAMILTONIAN_XC_FUNCTIONAL_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;	
	using namespace Catch::literals;
	using Catch::Approx;

	SECTION("LDA"){
		inq::hamiltonian::xc_functional ldafunctional(XC_LDA_X, 1);
		CHECK(ldafunctional.exx_coefficient() == 0.0);
	}

	SECTION("LSDA"){
		inq::hamiltonian::xc_functional ldafunctional(XC_LDA_X, 2);
		CHECK(ldafunctional.exx_coefficient() == 0.0);
	}
		
	SECTION("GGA"){
		inq::hamiltonian::xc_functional ggafunctional(XC_GGA_X_PBE, 1);
		CHECK(ggafunctional.exx_coefficient() == 0.0);
	}

	SECTION("Spin GGA"){
		inq::hamiltonian::xc_functional ggafunctional(XC_GGA_X_PBE, 2);
		CHECK(ggafunctional.exx_coefficient() == 0.0);
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

