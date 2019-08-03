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

#endif

