/* -*- indent-tabs-mode: t -*- */

#ifndef IONIC__SPECIES
#define IONIC__SPECIES

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <pseudopod/element.hpp>
#include <utils/merge_optional.hpp>
#include <vector>
#include <cmath>
#include <optional>

namespace inq {
namespace input {

class species : public pseudo::element {

	std::optional<std::string> symbol_;
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
		
	auto pseudo(const std::string & pseudo_file){
		species rspec = *this;
		rspec.pseudo_file_ = pseudo_file;
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

	template<class OStream>
	friend OStream & operator<<(OStream & out, species const & self){
		out << self.symbol();
		return out;
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

	SECTION("Constructor"){
		
		auto s = input::species("Xe");
		
		CHECK(s.atomic_number() == 54);
		CHECK(not s.has_file());
	}

	SECTION("Constructor with options"){
		
		auto s = input::species("Xe").mass(20);
		
		CHECK(s.atomic_number() == 54);
		CHECK(not s.has_file());
		CHECK(s.mass() == 36457.77_a);
	}

	SECTION("Option mass"){
		
		auto s = input::species("U").mass(235);
		
		CHECK(s.symbol() == "U");
		CHECK(s.mass() == 428378.7975_a);
	}
	
	SECTION("Option symbol"){
		
		auto s = input::species("U").symbol("U235").mass(235);
		
		CHECK(s.symbol() == "U235");
		CHECK(s.mass() == 428378.7975_a);
	}

	SECTION("Option pseudopotential"){
		
		auto s = input::species("He").pseudo("hola");
		
		CHECK(s.symbol() == "He");
		CHECK(s.has_file());
		CHECK(s.file_path() == "hola");
	}
	
	
}
#endif
