/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef SYSTEMS__IONS
#define SYSTEMS__IONS

#include <cfloat>

#include <ions/geometry.hpp>
#include <ions/unitcell.hpp>

namespace systems {

  class ions {

  public:

    template<class lattice_vectors_type>
    ions(const lattice_vectors_type & arg_lattice_vectors, const ::ions::geometry & geo_arg):
      cell_(arg_lattice_vectors),
      geo_(geo_arg){
      
      geo_.info(std::cout);
      cell_.info(std::cout);

    }

    auto & geo() const {
      return geo_;
    }

    auto & cell() const {
      return cell_;
    }
    
  private:
    
    ::ions::UnitCell cell_;
    ::ions::geometry geo_;

  };  
  
}

#endif

