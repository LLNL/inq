#ifndef PARSER_INPUTFILE
#define PARSER_INPUTFILE

#include <string>
#include <iostream>
#include <cstdlib>

#include "liboct_parser.h"

namespace parser {

  class input_file {
    
  public:

    input_file(const std::string & file_name){
      int dont_write = 0;
      int set_used = 1;      
      parse_init("parser.log", &dont_write);
      parse_input(file_name.c_str(), set_used);
      
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
      type dummy_default;
      if(defines(variable_name)){
	return parse(variable_name, dummy_default);
      } else {
	std::cerr << "Error: required variable " << variable_name << " is not defined in the input file" << std::endl;
	exit(1);
      }
    }

  private:
    
  };

}

#endif

// Local Variables:
// mode: c++
// coding: utf-8
// End:
