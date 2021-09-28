/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__MATH__COMPLEX
#define INQ__MATH__COMPLEX

#include <inq_config.h>

#include <complex>
#include <gpu/run.hpp>

#ifdef ENABLE_CUDA
#include <thrust/complex.h>
#endif

namespace inq {

#ifdef ENABLE_CUDA
using complex = thrust::complex<double>;
#else
using complex = std::complex<double>;
#endif

template <class Type>
struct compat_type {
	using type = Type;
};

template <>
struct compat_type<complex> {
	using type = std::complex<double>;
};

GPU_FUNCTION inline double real(const double & x){
	return x;
}

GPU_FUNCTION inline double real(const complex & z){
	return z.real();
}

GPU_FUNCTION inline double imag(const double &){
	return 0.0;
}

GPU_FUNCTION inline auto imag(const complex & z){
	return z.imag();
}

GPU_FUNCTION inline auto norm(const double & x){
	return x*x;
}

GPU_FUNCTION inline double conj(const double & x){
	return x;
}

GPU_FUNCTION inline auto mul(const double x, const double y){
	return x*y;
}

GPU_FUNCTION inline auto mul(const double x, const complex & z){
	return complex{x*real(z), x*imag(z)};
}

GPU_FUNCTION inline auto mul(const complex & z, const double x){
	return complex{real(z)*x, imag(z)*x};
}

GPU_FUNCTION inline auto mul(const complex & z1, const complex z2){
	return complex{real(z1)*real(z2) - imag(z1)*imag(z2), real(z1)*imag(z2) + imag(z1)*real(z2)};
}

GPU_FUNCTION inline auto fabs(complex const & z){
	return abs(z);
}

}

#ifdef INQ_MATH_COMPLEX_UNIT_TEST
#undef INQ_MATH_COMPLEX_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("Class math::complex", "[math::complex]"){
  
	using namespace inq;
	using namespace Catch::literals;

	static_assert(std::is_same<compat_type<double>::type, double>::value, "check blas types match");
	static_assert(std::is_same<compat_type<inq::complex>::type, std::complex<double>>::value, "check blas types match");

}

#endif
#endif


