/* -*- indent-tabs-mode: t -*- */

#ifndef MATH_COMPLEX
#define MATH_COMPLEX
#include <complex>

#ifdef HAVE_CUDA
#define GPU_FUNCTION __host__ __device__
#else
#define GPU_FUNCTION
#endif

using complex = std::complex<double>;

GPU_FUNCTION inline double conj(const double & x){
	return x;
}

GPU_FUNCTION inline complex conj(const complex & z){
	return complex(std::real(z), -std::imag(z));
}

GPU_FUNCTION inline double real(const double & x){
	return x;
}

GPU_FUNCTION inline double real(const complex & z){
	return std::real(z);
}

GPU_FUNCTION inline double imag(const double &){
	return 0.0;
}

GPU_FUNCTION inline auto imag(const complex & z){
	return std::imag(z);
}

GPU_FUNCTION inline auto operator*(const double x, const complex & z){
	return complex{x*real(z), x*imag(z)};
}

GPU_FUNCTION inline auto operator*(const complex & z, const double x){
	return complex{x*real(z), x*imag(z)};
}

GPU_FUNCTION inline auto operator+(const complex & z1, const complex z2){
	return complex{real(z1) + real(z2), imag(z1) + imag(z2)};
}

GPU_FUNCTION inline auto & operator+=(complex & z1, const complex z2){
	z1 = z1 + z2;
	return z1;
}

GPU_FUNCTION inline auto operator*(const complex & z1, const complex z2){
	return complex{real(z1)*real(z2) - imag(z1)*imag(z2), real(z1)*imag(z2) + imag(z1)*real(z2)};
}

GPU_FUNCTION inline auto norm(const double & x){
	return x*x;
}

inline auto norm(const complex & z){
	return std::norm(z);
}

#endif

