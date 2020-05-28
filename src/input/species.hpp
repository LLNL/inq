/* -*- indent-tabs-mode: t -*- */

#ifndef INPUT__SPECIES
#define INPUT__SPECIES

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

#include <pseudopod/element.hpp>
#include <utils/merge_optional.hpp>
#include <vector>
#include <cmath>
#include <nonstd/optional.hpp>

namespace inq {
namespace input {

class species : public pseudo::element {

	struct options;
		
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
		return opts.mass_.value_or(element::mass());
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

	struct options {
	private:
			
		nonstd::optional<std::string> symbol_;
		nonstd::optional<std::string> pseudo_file_;		
		nonstd::optional<double> mass_;
		nonstd::optional<bool> filter_;
			
		options(){
		}

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

#ifdef INQ_UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class input::species", "[input::species]") {
  
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
		CHECK(s.mass() == 20.0_a);
		
	}

	SECTION("Option mass"){
		
		input::species s = pseudo::element("U") | input::species::mass(235);
		
		CHECK(s.symbol() == "U");
		CHECK(s.mass() == 235.0_a);
		
	}
	
	SECTION("Option symbol"){
		
		input::species s = "U" | input::species::symbol("U235") | input::species::mass(235);
		
		CHECK(s.symbol() == "U235");
		CHECK(s.mass() == 235.0_a);
		
	}

	SECTION("Option pseudopotential"){
		
		input::species s = "He" | input::species::pseudo("hola");
		
		CHECK(s.symbol() == "He");
		CHECK(s.has_file());
		CHECK(s.file_path() == "hola");
		
	}
	
	
}


#endif


#endif
