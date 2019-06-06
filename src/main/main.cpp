#include <ions/geometry.hpp>
#include <ions/unitcell.hpp>
#include <basis/plane_wave.hpp>
#include <ions/geometry.hpp>
#include <hamiltonian/atomic_potential.hpp>

#include <multi/array.hpp>

#include <complex>

int main(){

  using math::d3vector;
  
  double ecut = 30.0;
  double ll = 10.0;
    
  ions::geometry geo(config::path::unit_tests_data() + "benzene.xyz");
  hamiltonian::atomic_potential pot(geo.number_of_atoms(), geo.atoms());

  ions::UnitCell cell(d3vector(ll, 0.0, 0.0), d3vector(0.0, ll, 0.0), d3vector(0.0, 0.0, ll));
  basis::plane_wave pw(cell, ecut);

  boost::multi::array<double, 3> density({pw.rsize()[0], pw.rsize()[1], pw.rsize()[2]});

  for(int ix = 0; ix < pw.rsize()[0]; ix++){
    for(int iy = 0; iy < pw.rsize()[1]; iy++){
      for(int iz = 0; iz < pw.rsize()[2]; iz++){
	density[ix][iy][iz] = 0.0;
      }
    }
  }

  density[0][0][0] = -1.0;
   
  
}
