/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__SOLVERS__FIRE
#define INQ__SOLVERS__FIRE

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <operations/sum.hpp>
#include <solvers/velocity_verlet.hpp>

#include <cassert>

namespace inq {
namespace solvers {

template <class ArrayType, class ForceFunction>
void fire(ArrayType & xx, double step, double tolforce, ForceFunction const & func){

  auto file = std::ofstream("fire.dat");
  
  auto alpha_start = 0.1;
  auto dt = step;
  auto alpha = alpha_start;
  auto p_times = 0;
  auto f_alpha = 0.99;
  auto n_min = 5;
  auto f_inc = 1.1;
  auto dt_max = 10.0*dt;
  auto f_dec = 0.5;
  auto const mass = 1.0;
  auto const maxiter = 200;

  auto vel = ArrayType(xx.size(), {0.0, 0.0, 0.0});
  for(int iiter = 0; iiter < maxiter; iiter++){

    auto force = func(xx);
    auto p_value = operations::sum(force, vel, [](auto fo, auto ve) { return dot(fo, ve);});

    std::cout << iiter << '\t' << xx[0][0] << '\t' << vel[0][0] << '\t' << force[0][0] << '\t' << p_value << '\t' << std::endl;
    file << iiter << '\t' << xx[0][0] << '\t' << vel[0][0] << '\t' << force[0][0] << '\t' << p_value << '\t' << dt << std::endl;
    
    auto max_force = 0.0;
    for(auto ii = 0; ii < force.size(); ii++) max_force = std::max(max_force, fabs(force[ii]));

    if(max_force < tolforce) break;
    
    auto norm_vel = operations::sum(vel, [](auto xx) { return norm(xx); });
    auto norm_force = operations::sum(force, [](auto xx) { return norm(xx); });
    for(auto ii = 0; ii < vel.size(); ii++) vel[ii] = (1.0 - alpha)*vel[ii] + alpha*force[ii]*sqrt(norm_vel/norm_force);
      
    if(p_value > 0.0) {
      if(p_times > n_min) {
        dt = std::min(dt*f_inc, dt_max);
        alpha *= f_alpha;
      }

      p_times++;
    } else {
      
      p_times = 0;
      dt *= f_dec;
      alpha = alpha_start;
      for(auto ii = 0; ii < vel.size(); ii++) vel[ii] = vector3{0.0, 0.0, 0.0};
    }

    for(auto ii = 0; ii < vel.size(); ii++) {
      vel[ii] += force[ii]*dt/mass;
      xx[ii]  += vel[ii]*dt;
    }
    
  }

}

}
}
#endif

///////////////////////////////////////////////////////////////////

#ifdef INQ_SOLVERS_FIRE_UNIT_TEST
#undef INQ_SOLVERS_FIRE_UNIT_TEST

#include <catch2/catch_all.hpp>

#include <gpu/array.hpp>

using namespace inq;
using namespace Catch::literals;
using Catch::Approx;

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {


  SECTION("Lennard Jones 2 atoms"){
    
    auto pos = gpu::array<vector3<double>, 1>(2);
    pos[0] = vector3{-1.0, 0.0, 0.0};
    pos[1] = vector3{ 1.0, 0.0, 0.0};

    auto lj2 = [] (auto xx) {
      auto dxx = xx[0] - xx[1];
      auto rr2 = norm(dxx);
      auto rr6 = rr2*rr2*rr2;
      auto rr12 = rr6*rr6;

      //      std::cout << 4.0*(1.0/rr12 - 1.0/rr6) << std::endl;
      
      auto force = gpu::array<vector3<double>, 1>(2);
      force[0] = 24.0/rr2*(2.0/rr12 - 1.0/rr6)*dxx;
      force[1] = -force[0];

      return force;
    };

    auto force = lj2(pos);

    CHECK(force[0][0] == Approx(0.181640625));
    CHECK(force[0][1] == Approx(0.0).margin(1e-12));
    CHECK(force[0][2] == Approx(0.0).margin(1e-12));
    CHECK(force[1][0] == Approx(-0.181640625));
    CHECK(force[1][1] == Approx(0.0).margin(1e-12));
    CHECK(force[1][2] == Approx(0.0).margin(1e-12));
  
    solvers::fire(pos, 0.1, 1e-5, lj2);
    
    CHECK(pos[0][0] == Approx(-0.5*pow(2, 1.0/6.0)));
    CHECK(pos[0][1] == Approx(0.0).margin(1e-12));
    CHECK(pos[0][2] == Approx(0.0).margin(1e-12));
    CHECK(pos[1][0] == Approx( 0.5*pow(2, 1.0/6.0)));
    CHECK(pos[1][1] == Approx(0.0).margin(1e-12));
    CHECK(pos[1][2] == Approx(0.0).margin(1e-12));

  }
	
}
#endif
