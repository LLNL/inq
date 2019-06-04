#include <ions/geometry.hpp>
#include <ions/unitcell.hpp>
#include <basis/basis.hpp>
#include <array.hpp>

#include <complex>

int main(){

  using math::d3vector;
  
  double ecut = 30.0;
  double ll = 10.0;

  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  basis::basis bas(cell, ecut);

  boost::multi::array<double, 3> density({bas.rsize()[0], bas.rsize()[1], bas.rsize()[2]});

  for(int ix = 0; ix < bas.rsize()[0]; ix++){
    for(int iy = 0; iy < bas.rsize()[1]; iy++){
      for(int iz = 0; iz < bas.rsize()[2]; iz++){
	density[ix][iy][iz] = 0.0;
      }
    }
  }

  density[0][0][0] = -1.0;
   
  
}
