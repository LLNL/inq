#ifndef IONS_POSITION
#define IONS_POSITION

#include <iostream>

#include <array>
#include <cmath>
#include <limits>

/*
  
  This class represents a position in space.

 */

namespace ions {

  class Position {

  public:

    Position(){}

    // constructs a Position from spherical coordinates
    static Position spherical(const double & r, const double & theta, const double & phi) {
      return Position(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));
    }
    
    Position(const double & x, const double & y, const double & z){
      coords[0] = x;
      coords[1] = y;
      coords[2] = z;
    }
    
    const Position & operator*=(const double & a){
      coords[0] *= a;
      coords[1] *= a;
      coords[2] *= a;
      return *this;
    }

    const Position & operator/=(const double & a){
      coords[0] /= a;
      coords[1] /= a;
      coords[2] /= a;
      return *this;
    }

    friend Position operator+(const Position & a, const Position & b){
      return Position(a.coords[0] + b.coords[0], a.coords[1] + b.coords[1], a.coords[2] + b.coords[2]);
    }
    
    friend Position operator-(const Position & a, const Position & b){
      return Position(a.coords[0] - b.coords[0], a.coords[1] - b.coords[1], a.coords[2] - b.coords[2]);
    }

    friend Position operator*(const double & a, const Position & pos){
      return Position(a*pos.coords[0], a*pos.coords[1], a*pos.coords[2]);
    }

    friend Position operator*(const Position & pos, const double & a){
      return a*pos;
    }

    friend Position operator/(const Position & pos, const double & a){
      return Position(pos.coords[0]/a, pos.coords[1]/a, pos.coords[2]/a);
    }
    
    friend std::istream & operator>>(std::istream & input, Position & pos){
      input >> pos.coords[0] >> pos.coords[1] >> pos.coords[2];
      return input;
    }

    friend std::ostream & operator<<(std::ostream & output, const Position & pos){
      output << ' ' << pos.coords[0] << ' ' << pos.coords[1] << ' ' << pos.coords[2];
      return output;
    }

    friend double norm_sq(const Position & pos){
      return pos.x()*pos.x() + pos.y()*pos.y() + pos.z()*pos.z();
    }    
    
    friend double norm(const Position & pos){
      return std::sqrt(norm_sq(pos));
    }

    double costheta() const {
      if(fabs(z()) < std::numeric_limits<double>::epsilon()) return 0.0;
      return z()/std::sqrt(norm_sq(*this));
    }

    double phi() const {
      return std::atan2(y(), x());
    }

    const double & x() const {
      return coords[0];
    }

    const double & y() const {
      return coords[1];
    }

    const double & z() const {
      return coords[2];
    }

  private:

    std::array<double, 3> coords;
    
  };

}

#endif

// Local Variables:
// mode: c++
// coding: utf-8
// End:
