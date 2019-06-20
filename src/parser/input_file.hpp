#ifndef PARSER_INPUTFILE
#define PARSER_INPUTFILE

#include <string>
#include <iostream>
#include <cstdlib>

#include "liboct_parser.h"

namespace parser {

  class input_file {
    
  public:

    enum class error {
      FILE_NOT_FOUND,
      VARIABLE_NOT_FOUND
    };
    
    input_file(const std::string & file_name){
      int dont_write = 0;
      int set_used = 1;      
      parse_init("parser.log", &dont_write);
      if(parse_input(file_name.c_str(), set_used)){
	throw error::FILE_NOT_FOUND;
      }
    }

    ~input_file(){
      parse_end();
    }

    bool defines(const std::string & variable_name) const {
      return parse_isdef(variable_name.c_str());
    }
    
    int parse(const std::string & variable_name, const int & default_value) const {
      return parse_int(variable_name.c_str(), default_value);
    }

    double parse(const std::string & variable_name, const double & default_value) const {
      return parse_double(variable_name.c_str(), default_value);
    }

    std::string parse(const std::string & variable_name, const std::string & default_value) const {
      return parse_string(variable_name.c_str(), default_value.c_str());
    }

    bool parse(const std::string & variable_name, const bool & default_value) const {
      // call the integer version of parse
      return parse(variable_name, (default_value?1:0)) == 1;
    }


    // this version doest not receive a default, but fails when the variable is not found
    template <class type>
    type parse(const std::string & variable_name) const {
      if(defines(variable_name)){
	return parse(variable_name, type());
      } else {
	throw error::VARIABLE_NOT_FOUND;
      }
    }
    
  private:
    
  };

}

#ifdef UNIT_TEST
#include <catch2/catch.hpp>
#include <cmath>
#include <config/path.hpp>

TEST_CASE("Class parser::input_file", "[input_file]") {

  using namespace Catch::literals;

  SECTION("Non-existing file"){
    REQUIRE_THROWS(parser::input_file("/a_file_that_doesnt_exists"));
  }

  SECTION("Parse file"){
    parser::input_file input(config::path::unit_tests_data() + "input_file");

    // read with default value, this value is returned if the variable doesn't exist
    REQUIRE(input.parse("thisvariabledoesntexist", 10) == 10);

    REQUIRE(input.parse("doublevar", -1.0) == 199.33_a);

    REQUIRE(input.parse<std::string>("stringvar") == "hola");

    REQUIRE_THROWS(input.parse<float>("anothetvariablethatdoesntexist"));

    REQUIRE(input.defines("intvar"));

    REQUIRE(input.parse<int>("intvar") == 100774);
 
  }

  
}
  
#endif

#endif

// Local Variables:
// eval:(setq indent-tabs-mode: t tab-width: 2)
// mode: c++
// coding: utf-8
// End:
