/* -*- indent-tabs-mode: t -*- */

#ifndef INPUT__SPECIES
#define INPUT__SPECIES

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

	class options;
		
public:

	species(const pseudo::element & arg_el, const options & arg_options = {}):
		pseudo::element(arg_el),
		opts(arg_options){
	}

	species(const species &) = default;	
	species(species &&) = default;
	species & operator=(const species &) = delete;	
	species & operator=(species &&) = delete;
		
	friend species operator|(const species & spec, const options & opts){
		auto rspec = spec;
		rspec.opts = rspec.opts | opts;
		return rspec;
	}

	friend species operator|(const std::string & arg_symbol, const options & opts){
		auto rspec = species(arg_symbol);
		rspec.opts = rspec.opts | opts;
		return rspec;
	}
		
	static auto symbol(const std::string & arg_symbol){
		options ropt;
		ropt.symbol_ = arg_symbol;
		return ropt;
	}
		
	static auto pseudo(const std::string & pseudo_file){
		options ropt;
		ropt.pseudo_file_ = pseudo_file;
		return ropt;
	}

	static auto mass(const double arg_mass){
		options ropt;
		ropt.mass_ = arg_mass;
		return ropt;
	}

	auto has_file() const {
		return opts.pseudo_file_.has_value();
	}

	auto const & file_path() const {
		return opts.pseudo_file_.value();
	}

	auto symbol() const {
		using pseudo::element;
		return opts.symbol_.value_or(element::symbol());
	}

	auto mass() const {
		using pseudo::element;
		return 1822.8885*opts.mass_.value_or(element::mass());
	}

	static auto nofilter() {
		options ropt;
		ropt.filter_ = false;
		return ropt;
	}

	static auto filter() {
		options ropt;
		ropt.filter_ = true;
		return ropt;
	}
		
	auto filter_pseudo() const {
		return opts.filter_.value_or(true);
	}
				
private:

	class options {
			
		std::optional<std::string> symbol_;
		std::optional<std::string> pseudo_file_;		
		std::optional<double> mass_;
		std::optional<bool> filter_;
			
	public:

		friend options operator|(const options & opt1, const options & opt2){
			options ropt;

			using inq::utils::merge_optional;
				
			ropt.symbol_ = merge_optional(opt1.symbol_, opt2.symbol_);
			ropt.pseudo_file_ = merge_optional(opt1.pseudo_file_, opt2.pseudo_file_);
			ropt.mass_ = merge_optional(opt1.mass_, opt2.mass_);
			ropt.filter_ = merge_optional(opt1.filter_, opt2.filter_);
				
			return ropt;
		}

		friend class species;
	};

	options opts;
		
};

}
}
#endif

#ifdef INQ_INPUT_SPECIES_UNIT_TEST
#undef INQ_INPUT_SPECIES_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

	SECTION("Constructor"){
		
		input::species s(pseudo::element("Xe"));
		
		CHECK(s.atomic_number() == 54);
		CHECK(not s.has_file());
	}

	SECTION("Constructor with options"){
		
		input::species s(pseudo::element("Xe"), input::species::mass(20));
		
		CHECK(s.atomic_number() == 54);
		CHECK(not s.has_file());
		CHECK(s.mass() == 36457.77_a);
	}

	SECTION("Option mass"){
		
		input::species s = pseudo::element("U") | input::species::mass(235);
		
		CHECK(s.symbol() == "U");
		CHECK(s.mass() == 428378.7975_a);
	}
	
	SECTION("Option symbol"){
		
		input::species s = "U" | input::species::symbol("U235") | input::species::mass(235);
		
		CHECK(s.symbol() == "U235");
		CHECK(s.mass() == 428378.7975_a);
	}

	SECTION("Option pseudopotential"){
		
		input::species s = "He" | input::species::pseudo("hola");
		
		CHECK(s.symbol() == "He");
		CHECK(s.has_file());
		CHECK(s.file_path() == "hola");
	}
	
	
}
#endif
