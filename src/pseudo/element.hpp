#ifndef PSEUDO_ELEMENT_HPP
#define PSEUDO_ELEMENT_HPP

/*
 Copyright (C) 2018 Xavier Andrade

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

#include <string>
#include <cassert>
#include <cctype>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <map>
#include <cstdlib>

#include <config/path.hpp>

namespace pseudo {
  
  class element {

  private:

    struct properties;
    typedef std::map<std::string, properties> map_type;
    
  public:
    
    element(const std::string & symbol = "none"):symbol_(trim(symbol)){
      symbol_[0] = std::toupper(symbol_[0]);
      for(unsigned ii = 1; ii < symbol_.size(); ii++) symbol_[ii] = std::tolower(symbol_[ii]);

      map(); //make sure the map is initialized
    }

    element(int atomic_number){

      //special case: avoid ambiguity between isotopes
      if(atomic_number == 1){
	symbol_ = 'H';
	return;
      }
      
      for(map_type::iterator it = map().begin(); it != map().end(); ++it){
	if(it->second.z_ == atomic_number) {
	  symbol_ = trim(it->first);
	  break;
	}
      }
    }
    
    bool valid() const {
      return map().find(symbol_) != map().end();
    }

    double charge() const {
      return -1.0*atomic_number();
    }

    const std::string & symbol() const{
      return symbol_;
    }
    
    int atomic_number() const {
      return map().at(symbol_).z_;
    }
    
    double mass() const{
      return map().at(symbol_).mass_;
    }

    double vdw_radius() const{
      return map().at(symbol_).vdw_radius_;
    }

    bool operator==(const element & el) const {
      return atomic_number() == el.atomic_number();
    }

    struct hash {
      auto operator()(const element & el){
	std::hash<std::string> str_hash;
	return str_hash(el.symbol());
      }
    };
    
  private:

    struct properties {
      int z_;
      double mass_;
      double vdw_radius_;
    };
    
    static map_type & map(){
      
      static map_type map;
      if(map.empty()){

	std::string filename = config::path::share() + "/pseudopotentials/elements.dat";
	
	std::ifstream file(filename.c_str());

	if(!file){
	  std::cerr << "Internal error: cannot open file '" << filename << "'." << std::endl;
	  exit(EXIT_FAILURE);
	}
	
	while(true){
	  std::string symbol;
	  
	  file >> symbol;

	  if(file.eof()) break;
	  
	  if(symbol[0] == '#'){
	    getline(file, symbol);
	    continue;
	  }
	  
	  properties prop;
	  
	  file >> prop.z_ >> prop.mass_ >> prop.vdw_radius_;

	  if(file.eof()) break;
	  
	  map[symbol] = prop;
	  
	}
      }

      return map;
    }
    
    static std::string & ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r "){
      str.erase(0, str.find_first_not_of(chars));
      return str;
    }
 
    static std::string & rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r "){
      str.erase(str.find_last_not_of(chars) + 1);
      return str;
    }

  public:
    
    static std::string trim(std::string str, const std::string& chars = "\t\n\v\f\r "){
      return ltrim(rtrim(str, chars), chars);
    }
    
  private:
    
    std::string symbol_;
  
  };

}


#ifdef UNIT_TEST
#include <catch2/catch.hpp>

TEST_CASE("Class pseudo::element", "[element]") {

  using Catch::Matchers::WithinULP;

  SECTION("Define element by atomic number"){
    pseudo::element el(27);

    REQUIRE(el == pseudo::element("Co"));
    REQUIRE(el.valid());
    REQUIRE(el.atomic_number() == 27);
    REQUIRE(el.symbol() == "Co");
    REQUIRE_THAT(el.charge(), WithinULP(-27.0, 1));
    REQUIRE_THAT(el.mass(), WithinULP(58.9331944, 1));
    REQUIRE_THAT(el.vdw_radius(), WithinULP(4.535342, 1));
  }
  
  SECTION("Define element by symbol"){
    pseudo::element el("I");

    REQUIRE(el == pseudo::element(53));
    REQUIRE(el.valid());
    REQUIRE(el.atomic_number() == 53);
    REQUIRE(el.symbol() == "I");
    REQUIRE_THAT(el.charge(), WithinULP(-53.0, 1));
    REQUIRE_THAT(el.mass(), WithinULP(126.904473, 1));
    REQUIRE_THAT(el.vdw_radius(), WithinULP(3.855041, 1));
  }
  
  SECTION("Invalid element by atomic number"){
    pseudo::element el(200);

    REQUIRE(!el.valid());
  }

  SECTION("Invalid element by symbol"){
    pseudo::element el("X");

    REQUIRE(!el.valid());
  }
  
}
#endif

#endif

// Local Variables:
// eval:(setq indent-tabs-mode t tab-width 2)
// mode: c++
// coding: utf-8
// End:
