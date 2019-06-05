#ifndef PSEUDOPOTENTIAL_HPP
#define PSEUDOPOTENTIAL_HPP

#include <pseudo/psml.hpp>
#include <pseudo/qso.hpp>
#include <pseudo/upf1.hpp>
#include <pseudo/upf2.hpp>
#include <pseudo/psp8.hpp>
#include <pseudo/detect_format.hpp>

namespace pseudo {

  class pseudopotential {

  public:

    enum class error {
      FILE_NOT_FOUND,
      UNKNOWN_FORMAT,
      UNSUPPORTED_FORMAT
    };
    
    pseudopotential(const std::string & filename){

      pseudo::format format = pseudo::detect_format(filename);
      
      if(format == pseudo::format::FILE_NOT_FOUND) throw error::FILE_NOT_FOUND;
      if(format == pseudo::format::UNKNOWN) throw error::UNKNOWN_FORMAT;
      
      std::cout << "  <!-- SpeciesReader opening file " << filename << " -->" << std::endl;
      
      switch(format){
      case pseudo::format::QSO:
	std::cout << "  <!--   format: QSO -->" << std::endl;
	pseudo_ = new pseudo::qso(filename);
	break;
      case pseudo::format::UPF1:
	std::cout << "  <!--   format: UPF1 -->" << std::endl;
	pseudo_ = new pseudo::upf1(filename, /*uniform_grid = */ true);
	break;
      case pseudo::format::UPF2:
	std::cout << "  <!--   format: UPF2 -->" << std::endl;
	pseudo_ = new pseudo::upf2(filename, /*uniform_grid = */ true);
	break;
      case pseudo::format::PSML:
	std::cout << "  <!--   format: PSML -->" << std::endl;
	pseudo_ = new pseudo::psml(filename, /*uniform_grid = */ true);
	break;
      case pseudo::format::PSP8:
	std::cout << "  <!--   format: PSP8 -->" << std::endl;
	pseudo_ = new pseudo::psp8(filename);
	break;
      default:
	throw error::UNSUPPORTED_FORMAT;
      }
      
      std::cout << "  <!--   size:   " << pseudo_->size() << " -->" << std::endl;
      
    }

    ~pseudopotential(){
      delete pseudo_;
    }
    
  private:
    
    base * pseudo_;    
    
  };
  
}

#ifdef UNIT_TEST
#include <catch.hpp>

TEST_CASE("class pseudo::pseudopotential", "[pseudopotential]") {
  
  using namespace Catch::literals;

  SECTION("Non-existing file"){
    REQUIRE_THROWS(pseudo::pseudopotential("/this_file_doesnt_exists"));
  }

  SECTION("Non-pseudopotential file"){
    REQUIRE_THROWS(pseudo::pseudopotential(SHARE_DIR + std::string("/unit_tests_data/benzene.xyz")));
  }

  SECTION("Qbox pseudopotential file"){
    pseudo::pseudopotential ps(SHARE_DIR + std::string("/unit_tests_data/I_HSCV_LDA-1.0.xml"));
  }

  SECTION("UPF2 pseudopotential file"){
    pseudo::pseudopotential ps(SHARE_DIR + std::string("/unit_tests_data/W_ONCV_PBE-1.0.upf"));
  }

  SECTION("UPF1 pseudopotential file"){
    pseudo::pseudopotential ps(SHARE_DIR + std::string("/unit_tests_data/F.UPF"));
  }

  SECTION("UPF2 pseudopotential file"){
    pseudo::pseudopotential ps(SHARE_DIR + std::string("/unit_tests_data/78_Pt_r.oncvpsp.psp8"));
  }

}
#endif
 
#endif
