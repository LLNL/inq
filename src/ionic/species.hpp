/* -*- indent-tabs-mode: t -*- */

#ifndef IONIC__SPECIES
#define IONIC__SPECIES

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <pseudopod/element.hpp>
#include <pseudopod/set_id.hpp>
#include <utils/merge_optional.hpp>
#include <vector>
#include <cmath>
#include <optional>

namespace inq {
namespace ionic {

class species : public pseudo::element {

	std::optional<std::string> symbol_;
	std::optional<pseudo::set_id> pseudo_set_;
	std::optional<std::string> pseudo_file_;	
	std::optional<double> mass_;
	std::optional<bool> filter_;

public:
	
	species(char const * const arg_symbol):
		pseudo::element(std::string(arg_symbol)){
	}
	
	species(std::string const & arg_symbol):
		pseudo::element(arg_symbol){
	}

	species(int atomic_number):
		pseudo::element(atomic_number){
	}
	
	auto symbol(const std::string & arg_symbol) const{
		species rspec = *this;
		rspec.symbol_ = arg_symbol;
		return rspec;
	}
		
	auto pseudo_file(const std::string & file){
		species rspec = *this;
		rspec.pseudo_file_ = file;
		return rspec;
	}

	auto pseudo_set(pseudo::set_id const & set){
		species rspec = *this;
		rspec.pseudo_set_ = set;
		return rspec;
	}
	
	auto mass(const double arg_mass){
		species rspec = *this;
		rspec.mass_ = arg_mass;
		return rspec;
	}

	auto has_file() const {
		return pseudo_file_.has_value();
	}

	auto const & file_path() const {
		return pseudo_file_.value();
	}
	
	auto has_pseudo_set() const {
		return pseudo_set_.has_value();
	}

	auto const & pseudo_set() const {
		return pseudo_set_.value();
	}

	auto symbol() const {
		using pseudo::element;
		return symbol_.value_or(element::symbol());
	}

	auto mass() const {
		using pseudo::element;		
		return 1822.8885*mass_.value_or(element::mass());
	}

	auto nofilter() {
		species rspec = *this;
		rspec.filter_ = false;
		return rspec;
	}

	auto filter() {
		species rspec = *this;
		rspec.filter_ = true;
		return rspec;
	}
		
	auto filter_pseudo() const {
		return filter_.value_or(true);
	}

	void save(parallel::communicator & comm, std::string const & dirname) const {
		auto error_message = "INQ error: Cannot save species to directory '" + dirname + "'.";

		utils::create_directory(comm, dirname);
		
		utils::save_value   (comm, dirname + "/element_symbol", element::symbol(), error_message);
		utils::save_optional(comm, dirname + "/symbol",         symbol_,           error_message);
		utils::save_optional(comm, dirname + "/pseudo_set",     pseudo_set_,       error_message);
		utils::save_optional(comm, dirname + "/pseudo_file",    pseudo_file_,      error_message);
		utils::save_optional(comm, dirname + "/mass",           mass_,             error_message);
		utils::save_optional(comm, dirname + "/filter",         filter_,           error_message);
	}
	
	static auto load(std::string const & dirname) {

		auto error_message = "INQ error: Cannot load species from directory '" + dirname + "'.";

		std::string element_symbol;
		utils::load_value   (dirname + "/element_symbol", element_symbol, error_message);

		auto read_species = species{element_symbol};
		utils::load_optional(dirname + "/symbol",         read_species.symbol_     );
		utils::load_optional(dirname + "/pseudo_set",     read_species.pseudo_set_ );
		utils::load_optional(dirname + "/pseudo_file",    read_species.pseudo_file_);
		utils::load_optional(dirname + "/mass",           read_species.mass_       );
		utils::load_optional(dirname + "/filter",         read_species.filter_     );

		return read_species;
	}
	
};

}
}
#endif

#ifdef INQ_IONIC_SPECIES_UNIT_TEST
#undef INQ_IONIC_SPECIES_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};

	SECTION("Constructor"){
		
		auto s = ionic::species("Xe");
		
		CHECK(s.atomic_number() == 54);
		CHECK(not s.has_file());
		CHECK(not s.has_pseudo_set());
		
		s.save(comm, "save_species_xe");
		auto read_s = ionic::species::load("save_species_xe");

		CHECK(read_s.atomic_number() == 54);
		CHECK(not read_s.has_file());
		CHECK(not read_s.has_pseudo_set());
		
	}

	SECTION("Constructor with options"){
		
		auto s = ionic::species("Xe").mass(20);
		
		CHECK(s.atomic_number() == 54);
		CHECK(not s.has_file());
		CHECK(s.mass() == 36457.77_a);
		CHECK(not s.has_pseudo_set());

		s.save(comm, "save_species_xe_mass");
		auto read_s = ionic::species::load("save_species_xe_mass");

		CHECK(read_s.atomic_number() == 54);
		CHECK(not read_s.has_file());
		CHECK(read_s.mass() == 36457.77_a);
		CHECK(not read_s.has_pseudo_set());
		
	}

	SECTION("Option mass"){
		
		auto s = ionic::species("U").mass(235);
		
		CHECK(s.symbol() == "U");
		CHECK(s.mass() == 428378.7975_a);
		CHECK(not s.has_pseudo_set());

		s.save(comm, "save_species_u");
		auto read_s = ionic::species::load("save_species_u");

		CHECK(read_s.symbol() == "U");
		CHECK(read_s.mass() == 428378.7975_a);
		CHECK(not read_s.has_pseudo_set());
		
	}
	
	SECTION("Option symbol"){
		
		auto s = ionic::species("U").symbol("U235").mass(235);
		
		CHECK(s.symbol() == "U235");
		CHECK(s.mass() == 428378.7975_a);
		CHECK(not s.has_pseudo_set());
		
		s.save(comm, "save_species_u235");
		auto read_s = ionic::species::load("save_species_u235");

		CHECK(read_s.symbol() == "U235");
		CHECK(read_s.mass() == 428378.7975_a);
		CHECK(not read_s.has_pseudo_set());
		
	}

	SECTION("Option pseudopotential"){
		
		auto s = ionic::species("He").pseudo_file("hola");
		
		CHECK(s.symbol() == "He");
		CHECK(s.has_file());
		CHECK(s.file_path() == "hola");
		CHECK(not s.has_pseudo_set());

		s.save(comm, "save_species_he_hola");
		auto read_s = ionic::species::load("save_species_he_hola");

		CHECK(read_s.symbol() == "He");
		CHECK(read_s.has_file());
		CHECK(read_s.file_path() == "hola");
		CHECK(not read_s.has_pseudo_set());
		
	}
	
	SECTION("Option pseudopotential set"){
		
		auto s = ionic::species("He").pseudo_set(pseudo::set_id::ccecp());
		
		CHECK(s.symbol() == "He");
		CHECK(not s.has_file());
		CHECK(s.has_pseudo_set());
		CHECK(s.pseudo_set() == pseudo::set_id::ccecp());

		s.save(comm, "save_species_he_ccecp");
		auto read_s = ionic::species::load("save_species_he_ccecp");
		
		CHECK(read_s.symbol() == "He");
		CHECK(not read_s.has_file());
		CHECK(read_s.has_pseudo_set());
		CHECK(read_s.pseudo_set() == pseudo::set_id::ccecp());
		
	}
	
}
#endif
