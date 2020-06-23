/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATH__COMPLEX
#define INQ__MATH__COMPLEX

#ifdef HAVE_CONFIG_H
#include <inq_config.h>
#endif

#include <complex>
#include <gpu/run.hpp>

namespace inq {

/*
// This is currently disabled since it doesn't work with multi::blas

#ifdef HAVE_CUDA
#include <thrust/complex.h>

using complex = thrust::complex<double>;
*/

using complex = std::complex<double>;

//real

GPU_FUNCTION inline double real(const double & x){
	return x;
}

GPU_FUNCTION inline double real(const complex & z){
	return z.real();
}

//imag

GPU_FUNCTION inline double imag(const double &){
	return 0.0;
}


GPU_FUNCTION inline auto imag(const complex & z){
	return z.imag();
}

//norm

GPU_FUNCTION inline auto norm(const double & x){
	return x*x;
}

GPU_FUNCTION inline auto norm(const complex & z){
	return real(z)*real(z) + imag(z)*imag(z);
}

GPU_FUNCTION inline auto fabs(const complex & z){
	return hypot(real(z), imag(z));
}

//conj

GPU_FUNCTION inline double conj(const double & x){
	return x;
}

#ifdef HAVE_CUDA
GPU_FUNCTION inline complex conj(const complex & z){
	return complex(real(z), -imag(z));
}
#endif

#ifdef HAVE_CUDA
// sum
GPU_FUNCTION inline auto operator+(const complex & z1, const complex z2){
	return complex{real(z1) + real(z2), imag(z1) + imag(z2)};
}

GPU_FUNCTION inline auto operator-(const complex & z1, const complex z2){
	return complex{real(z1) - real(z2), imag(z1) - imag(z2)};
}
#endif

GPU_FUNCTION inline auto & operator+=(complex & z1, const complex z2){
	z1 = z1 + z2;
	return z1;
}

GPU_FUNCTION inline auto & operator-=(complex & z1, const complex z2){
	z1 = z1 - z2;
	return z1;
}

// multiplication

#ifdef HAVE_CUDA
GPU_FUNCTION inline auto operator*(const double x, const complex & z){
	return complex{x*real(z), x*imag(z)};
}

GPU_FUNCTION inline auto operator*(const complex & z, const double x){
	return complex{real(z)*x, imag(z)*x};
}

GPU_FUNCTION inline auto operator*(const complex & z1, const complex z2){
	return complex{real(z1)*real(z2) - imag(z1)*imag(z2), real(z1)*imag(z2) + imag(z1)*real(z2)};
}
#endif

//division

GPU_FUNCTION inline auto operator/(const complex & z, const int ii){
	return complex{real(z)/ii, imag(z)/ii};
}

#ifdef HAVE_CUDA

GPU_FUNCTION inline auto operator/(const complex & z, const double x){
	return complex{real(z)/x, imag(z)/x};
}

GPU_FUNCTION inline auto operator/(const complex & x, const complex & y){
  using std::abs;

  double s = abs(y.real()) + abs(y.imag());

  double oos = double(1.0) / s;

  double ars = x.real() * oos;
  double ais = x.imag() * oos;
  double brs = y.real() * oos;
  double bis = y.imag() * oos;

  s = (brs * brs) + (bis * bis);

  oos = double(1.0) / s;

  complex quot( ((ars * brs) + (ais * bis)) * oos, ((ais * brs) - (ars * bis)) * oos);
  return quot;

}
#endif

}
#endif

