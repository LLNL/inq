/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

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

#include <pseudo/element.hpp>
#include <vector>
#include <cmath>

namespace input {

	class species : public pseudo::element {

		class options;
		
  public:
		
		species(const pseudo::element & arg_el, const options & arg_options = {}):
			pseudo::element(arg_el),
			opts(arg_options){
		}

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
			return opts.symbol_.value_or(pseudo::element::symbol());
		}

		auto mass() const {
			return opts.mass_.value_or(pseudo::element::mass());
		}
		
	private:

		struct options {
		private:
			
			std::optional<std::string> symbol_;
			std::optional<std::string> pseudo_file_;		
			std::optional<double> mass_;
			
			options(){
			}
			
			template <class opt_type>
			static opt_type merge_option(const opt_type & option1, const opt_type & option2){
				if(option2) return option2;
				if(option1) return option1;
				return opt_type{};
			}
			
		public:

			friend options operator|(const options & opt1, const options & opt2){
				options ropt;
				
				ropt.symbol_ = merge_option(opt1.symbol_, opt2.symbol_);
				ropt.pseudo_file_ = merge_option(opt1.pseudo_file_, opt2.pseudo_file_);
				ropt.mass_ = merge_option(opt1.mass_, opt2.mass_);
				
				return ropt;
			}

			friend class species;
			
		};

		options opts;
		
  };
  
}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("class input::species", "[input::species]") {
  
  using namespace Catch::literals;

	SECTION("Constructor"){
		
		input::species s(pseudo::element("Xe"));
		
		REQUIRE(s.atomic_number() == 54);
		REQUIRE(not s.has_file());
		
	}

	SECTION("Option mass"){
		
		input::species s = pseudo::element("U") | input::species::mass(235);
		
		REQUIRE(s.symbol() == "U");
		REQUIRE(s.mass() == 235.0_a);
		
	}
	
	SECTION("Option symbol"){
		
		input::species s = "U" | input::species::symbol("U235") | input::species::mass(235);
		
		REQUIRE(s.symbol() == "U235");
		REQUIRE(s.mass() == 235.0_a);
		
	}

	SECTION("Option pseudopotential"){
		
		input::species s = "He" | input::species::pseudo("hola");
		
		REQUIRE(s.symbol() == "He");
		REQUIRE(s.has_file());
		REQUIRE(s.file_path() == "hola");
		
	}
	
	
}


#endif


#endif

// Local Variables:
// eval:(setq indent-tabs-mode: t tab-width: 2)
// mode: c++
// coding: utf-8
// End:
