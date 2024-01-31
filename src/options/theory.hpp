/* -*- indent-tabs-mode: t -*- */

#ifndef OPTIONS__THEORY
#define OPTIONS__THEORY

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <utils/load_save.hpp>

#include <cassert>
#include <optional>
#include <stdexcept>

#include <hamiltonian/xc_functional.hpp>

namespace inq {
namespace options {

class theory {

public:

	enum class exchange_functional {
		NONE = XC_NONE,
		LDA = XC_LDA_X,
		PBE = XC_GGA_X_PBE,
		RPBE = XC_GGA_X_RPBE,		
		B = XC_GGA_X_B88,
		B3LYP = XC_HYB_GGA_XC_B3LYP,
		PBE0 = XC_HYB_GGA_XC_PBEH,
		HARTREE_FOCK = XC_HARTREE_FOCK
	};

	enum class correlation_functional {
		NONE = XC_NONE,
		LDA_PZ = XC_LDA_C_PZ,
		PBE = XC_GGA_C_PBE,
		LYP = XC_GGA_C_LYP
	};
	
private:

	std::optional<bool> hartree_potential_;
	std::optional<exchange_functional> exchange_;
	std::optional<correlation_functional> correlation_;
	std::optional<double> alpha_;

public:
	
	auto non_interacting() const {
		theory inter = *this;
		inter.hartree_potential_ = false;
		inter.exchange_ = exchange_functional::NONE;
		inter.correlation_ = correlation_functional::NONE;		
		return inter;
	}

	auto dft() const {
		theory inter = *this;
		inter.hartree_potential_ = true;
		return inter;
	}
	
	auto lda() const {
		theory inter = *this;
		inter.hartree_potential_ = true;
		inter.exchange_ = exchange_functional::LDA;
		inter.correlation_ = correlation_functional::LDA_PZ;		
		return inter;
	}

	auto hartree() const {
		theory inter = *this;
		inter.hartree_potential_ = true;
		inter.exchange_ = exchange_functional::NONE;
		inter.correlation_ = correlation_functional::NONE;		
		return inter;
	}
	
	auto hartree_fock() const {
		theory inter = *this;
		inter.hartree_potential_ = true;
		inter.exchange_ = exchange_functional::HARTREE_FOCK;
		inter.correlation_ = correlation_functional::NONE;		
		return inter;
	}
		
	auto exchange() const {
		return exchange_.value_or(exchange_functional::PBE);
	}

	auto correlation() const {
		return correlation_.value_or(correlation_functional::PBE);
	}

	auto pbe() const {
		theory inter = *this;
		inter.hartree_potential_ = true;
		inter.exchange_ = exchange_functional::PBE;
		inter.correlation_ = correlation_functional::PBE;
		return inter;
	}

	auto rpbe() const {
		theory inter = *this;
		inter.hartree_potential_ = true;
		inter.exchange_ = exchange_functional::RPBE;
		inter.correlation_ = correlation_functional::PBE;
		return inter;
	}

	auto pbe0()  const {
		theory inter = *this;
		inter.hartree_potential_ = true;		
		inter.exchange_ = exchange_functional::PBE0;
		inter.correlation_ = correlation_functional::NONE;
		return inter;
	}

	auto b3lyp()  const {
		theory inter = *this;
		inter.hartree_potential_ = true;		
		inter.exchange_ = exchange_functional::B3LYP;
		inter.correlation_ = correlation_functional::NONE;
		return inter;
	}

	auto exchange_coefficient() const {
		if(exchange() == exchange_functional::HARTREE_FOCK) return 1.0;
		if(exchange() == exchange_functional::NONE) return 0.0;
		throw std::runtime_error("inq internal error: exchange coefficient not known here for true functionals");
	}

	auto hartree_potential() const {
		return hartree_potential_.value_or(true);
	}
	
	auto self_consistent() const {
		return hartree_potential() or exchange() != exchange_functional::NONE or correlation() != correlation_functional::NONE;
	}

	auto induced_vector_potential(const double alpha = -4.0*M_PI){
		theory inter = *this;
		inter.alpha_ = alpha;
		return inter;
	}
	
	auto has_induced_vector_potential() const {
		return alpha_.has_value();
	}

	auto alpha_value() const {
		assert(alpha_.has_value());
		return alpha_.value();
	}
	
	void save(parallel::communicator & comm, std::string const & dirname) const {
		auto error_message = "INQ error: Cannot save theory to directory '" + dirname + "'.";

		comm.barrier();

		auto exception_happened = true;
		if(comm.root()) {
			
			try { std::filesystem::create_directories(dirname); }
			catch(...) {
				comm.broadcast_value(exception_happened);
				throw std::runtime_error(error_message);
			}

			utils::save_optional(comm, dirname + "/hartree_potential", hartree_potential_, error_message);
			utils::save_optional_enum(comm, dirname + "/exchange", exchange_, error_message);
			utils::save_optional_enum(comm, dirname + "/correlation", correlation_, error_message);
			utils::save_optional(comm, dirname + "/alpha", alpha_, error_message);

			exception_happened = false;
			comm.broadcast_value(exception_happened);
			
		} else {
			comm.broadcast_value(exception_happened);
			if(exception_happened) throw std::runtime_error(error_message);
		}
		
		comm.barrier();
	}
		
	static auto load(std::string const & dirname) {
		theory opts;
		
		utils::load_optional(dirname + "/hartree_potential", opts.hartree_potential_);
		utils::load_optional_enum(dirname + "/exchange", opts.exchange_);
		utils::load_optional_enum(dirname + "/correlation", opts.correlation_);
		utils::load_optional(dirname + "/alpha", opts.alpha_);
		
		return opts;
	}


	template<class OStream>
	friend OStream& operator<<(OStream& out, theory const & self){
		out << "Theory:\n";

		if(not self.hartree_potential() and self.exchange() == exchange_functional::NONE and self.correlation() == correlation_functional::NONE){
			out << " Non-interacting electrons" << std::endl;
			return out;
		}

		if(self.hartree_potential() and self.exchange() == exchange_functional::NONE and self.correlation() == correlation_functional::NONE){
			out << " Hartree (with self-interaction)\n\n";
			out << " [1] D. R. Hartree, Math. Proc. Camb. Philos. Soc. 24 1, 111 (1928)" << std::endl;
			return out;
		}

		if(self.exchange() != exchange_functional::NONE) {
			auto e_func = hamiltonian::xc_functional(int(self.exchange()), 1);
			
			out << "  "   << e_func.kind_name() << ":\n";
			out << "    " << e_func.family_name() << " - " << e_func.name() << "\n\n";
			out << e_func.references("    ") << "\n";		
		}

		if(self.correlation() != correlation_functional::NONE) {
			auto c_func = hamiltonian::xc_functional(int(self.correlation()), 1);
			
			out << "  "   << c_func.kind_name() << ":\n";
			out << "    " << c_func.family_name() << " - " << c_func.name() << "\n\n";
			out << c_func.references("    ") << "\n";
		}
		
		return out;
	}
	
};
    
}
}
#endif

#ifdef INQ_OPTIONS_THEORY_UNIT_TEST
#undef INQ_OPTIONS_THEORY_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace Catch::literals;

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	
	SECTION("Defaults"){

    options::theory inter;

		CHECK(inter.hartree_potential() == true);
		CHECK(inter.exchange() == options::theory::exchange_functional::PBE);
		CHECK(inter.correlation() == options::theory::correlation_functional::PBE);
		CHECK_THROWS(inter.exchange_coefficient());
		
		inter.save(comm, "theory_save_default");
		auto read_inter = options::theory::load("theory_save_non_default");

		CHECK(read_inter.hartree_potential() == true);
		CHECK(read_inter.exchange() == options::theory::exchange_functional::PBE);
		CHECK(read_inter.correlation() == options::theory::correlation_functional::PBE);
		CHECK_THROWS(read_inter.exchange_coefficient());
	}

  SECTION("Non interacting"){

    auto inter = options::theory{}.non_interacting();

		CHECK(inter.hartree_potential() == false);
		CHECK(inter.self_consistent() == false);
		CHECK(inter.exchange_coefficient() == 0.0);
		CHECK(inter.has_induced_vector_potential() == false);

		inter.save(comm, "theory_save_non_interacting");
		auto read_inter = options::theory::load("theory_save_non_interacting");

		CHECK(read_inter.hartree_potential() == false);
		CHECK(read_inter.self_consistent() == false);
		CHECK(read_inter.exchange_coefficient() == 0.0);
		CHECK(read_inter.has_induced_vector_potential() == false);
  }
	
  SECTION("Hartee-Fock"){

    auto inter = options::theory{}.hartree_fock();
		CHECK(inter.exchange_coefficient() == 1.0);
		CHECK(inter.has_induced_vector_potential() == false);

		inter.save(comm, "theory_save_hartree_fock");
		auto read_inter = options::theory::load("theory_save_hartree_fock");
		
		CHECK(read_inter.exchange_coefficient() == 1.0);
		CHECK(read_inter.has_induced_vector_potential() == false);
		
  }

	SECTION("Induced vector potential Yabana"){
		auto inter = options::theory{}.induced_vector_potential();
		CHECK(inter.has_induced_vector_potential() == true);
		CHECK(inter.alpha_value() == -4.0*M_PI);
		
		inter.save(comm, "theory_save_yabana");
		auto read_inter = options::theory::load("theory_save_yabana");

		CHECK(read_inter.has_induced_vector_potential() == true);
		CHECK(read_inter.alpha_value() == -4.0*M_PI);
	}

	SECTION("Induced vector potential Ullrich"){
		auto inter = options::theory{}.induced_vector_potential(0.2);
		CHECK(inter.has_induced_vector_potential() == true);
		CHECK(inter.alpha_value() == 0.2);
		
		inter.save(comm, "theory_save_ullrich");
		auto read_inter = options::theory::load("theory_save_ullrich");

		CHECK(read_inter.has_induced_vector_potential() == true);
		CHECK(read_inter.alpha_value() == 0.2);
	}

}
#endif
