/* -*- indent-tabs-mode: t -*- */

#ifndef D3VECTOR_H
#define D3VECTOR_H

////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
//
// Produced at the Lawrence Livermore National Laboratory.
// Written by Xavier Andrade <xavier@llnl.gov>, Erik Draeger
// <draeger1@llnl.gov> and Francois Gygi <fgygi@ucdavis.edu>.
// Based on the Qbox code by Francois Gygi Copyright (c) 2008
// LLNL-CODE-635376. All rights reserved. 
//
// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// d3vector.h
//
// double 3-vectors
//
////////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cmath>
#include <cassert>

#include <utils/gpu.hpp>

namespace math {
  class d3vector {
  public:
    
    double x, y, z;
    
    // explicit constructor to avoid implicit conversion from double to d3vector
    GPU_FUNCTION d3vector(double xv, double yv, double zv) : x{xv}, y{yv}, z{zv}{}
    GPU_FUNCTION explicit d3vector() : x(0.0), y(0.0), z(0.0) {}

    GPU_FUNCTION explicit d3vector(const double & vv) : x(vv), y(vv), z(vv) {}    

    GPU_FUNCTION explicit d3vector(const double* r) : x(r[0]), y(r[1]), z(r[2]) {}

    GPU_FUNCTION double& operator[](int i){
			static_assert(sizeof(*this) == sizeof(double)*3, "must be compatible with double[3]");
			return reinterpret_cast<double*>(this)[i];
    }
		
		GPU_FUNCTION double const& operator[](int i) const{
			static_assert(sizeof(*this) == sizeof(double)*3, "must be compatible with double[3]");
			return reinterpret_cast<double const*>(this)[i];
    }
		
    GPU_FUNCTION bool operator==(const d3vector &rhs) const
    {
      return x == rhs.x && y == rhs.y && z == rhs.z;
    }

    GPU_FUNCTION bool operator!=(const d3vector &rhs) const
    {
      return x != rhs.x || y != rhs.y || z != rhs.z;
    }

    GPU_FUNCTION d3vector& operator += ( const d3vector& rhs )
    {
      x += rhs.x; y += rhs.y; z += rhs.z;
      return *this;
    }

    GPU_FUNCTION d3vector& operator -= ( const d3vector& rhs )
    {
      x -= rhs.x; y -= rhs.y; z -= rhs.z;
      return *this;
    }

    GPU_FUNCTION d3vector& operator *= ( const double& rhs )
    {
      x *= rhs; y *= rhs; z *= rhs;
      return *this;
    }

    GPU_FUNCTION d3vector& operator /= ( const double& rhs )
    {
      x /= rhs; y /= rhs; z /= rhs;
      return *this;
    }

    GPU_FUNCTION friend const d3vector operator + (const d3vector& lhs, const d3vector& rhs )
    {
      return d3vector(lhs) += rhs;
    }

    GPU_FUNCTION friend const d3vector operator - ( const d3vector& a, const d3vector& b )
    {
      return d3vector(a) -= b;
    }

    GPU_FUNCTION friend d3vector operator - ( const d3vector& a ) // unary minus
    {
      return d3vector( -a.x, -a.y, -a.z );
    }

    GPU_FUNCTION friend d3vector operator * ( const double& a, const d3vector& b )
    {
      return d3vector(b) *= a;
    }

    GPU_FUNCTION friend d3vector operator * ( const d3vector& a, const double& b )
    {
      return d3vector(a) *= b;
    }

    GPU_FUNCTION friend d3vector operator / ( const d3vector& a, const double& b )
    {
      return d3vector(a) /= b;
    }

    // scalar product
    GPU_FUNCTION friend double operator * ( const d3vector& a, const d3vector& b )
    {
      return a.x * b.x + a.y * b.y + a.z * b.z ;
    }

    GPU_FUNCTION friend d3vector operator ^ ( const d3vector& a, const d3vector& b )
    {
      return d3vector( a.y * b.z - a.z * b.y ,
		       a.z * b.x - a.x * b.z ,
		       a.x * b.y - a.y * b.x );
    }

    friend d3vector rotate ( const d3vector& x, const d3vector& w )
    {
      if ( length(x) == 0.0 ) return x; // x has zero length
      double theta = length( w );       // rotate by zero
      if ( theta == 0.0 ) return x;
      d3vector ew = normalized ( w );
      d3vector v = w ^ x;
      if ( length( v ) == 0.0 ) return x; // x is parallel to the rotation axis
      v = normalized( v );
      d3vector u = v ^ ew;
      double p = x * u;
      return  (x*ew)*ew + p*cos(theta)*u + p*sin(theta)*v ;
    }

    GPU_FUNCTION friend double length( const d3vector& a )
    {
      return sqrt( a.x * a.x + a.y * a.y + a.z * a.z );
    }

    GPU_FUNCTION friend double norm( const d3vector& a )
    {
      return a.x * a.x + a.y * a.y + a.z * a.z;
    }

    GPU_FUNCTION friend d3vector normalized( const d3vector a )
    {
      return a / length( a );
    }

    friend std::ostream& operator << ( std::ostream& os, const d3vector& v )
    {
      os << v.x << " " << v.y << " " << v.z;
      return os;
    }

    friend std::istream& operator >> ( std::istream& is, d3vector& v )
    {
      is >> v.x >> v.y >> v.z ;
      return is;
    }

  };
}
#endif

