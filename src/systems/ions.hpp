/* -*- indent-tabs-mode: t; tab-width: 2 -*- */

#ifndef SYSTEMS__IONS
#define SYSTEMS__IONS

#include <cfloat>

#include <ions/geometry.hpp>
#include <ions/unitcell.hpp>

namespace systems {

  class ions {

  public:
    
    ions(const ::ions::geometry & geo_arg, const ::ions::UnitCell & arg_cell):
      geo_(geo_arg),
      cell_(arg_cell){
    }

    auto geo(){
      return geo_;
    }

    auto cell(){
      return cell_;
    }
    
  private:
    
    ::ions::geometry geo_;
    ::ions::UnitCell cell_;

  };  
  
}

#endif

