/* -*- indent-tabs-mode: t -*- */

#ifndef MATH_COMPLEX
#define MATH_COMPLEX
#include <complex>

using complex = std::complex<double>;

inline double conj(const double & x){
	return x;
}

inline complex conj(const complex & z){
	return std::conj(z);
}

inline double real(const double & x){
	return x;
}

inline double real(const complex & z){
	return std::real(z);
}

inline double imag(const double &){
	return 0.0;
}

inline double imag(const complex & z){
	return std::imag(z);
}

#endif

