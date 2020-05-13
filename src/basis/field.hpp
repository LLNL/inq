#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
OMPI_CXX=$CXX ../../blds/gcc/scripts/inc++ -x c++ $0 -o $0x&&$0x&&rm $0x;exit
#endif

#ifndef BASIS_FIELD
#define BASIS_FIELD

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo A. Correa

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <math/array.hpp>
#include <tinyformat/tinyformat.h>
#include <algorithm>
#include <utils/skeleton_wrapper.hpp>
#include <basis/real_space.hpp>
#include <math/complex.hpp>

#include <mpi3/environment.hpp>
#include <multi/adaptors/blas.hpp>

#include <fstream>

namespace basis {
	
	template<class Basis, typename T>
	class field {

	public:

		using element_type = T;
		using basis_type = Basis;
		using internal_array_type = math::array<element_type, 1>;
		
		field(const basis_type & basis, boost::mpi3::communicator & comm = boost::mpi3::environment::get_self_instance()):
			basis_comm_(comm),
			linear_(basis.part().local_size()),
			basis_(basis){

			assert(basis_.part().comm_size() == basis_comm_.size());
		}

		template <class OtherT>
		field(skeleton_wrapper<field<basis_type, OtherT>> const & skeleton)
			:field(skeleton.base.basis(), skeleton.base.basis_comm()){
		}

		template<class, class> friend class field;
		template<typename OtherT>
		field(field<basis_type, OtherT> const& o) 
		: basis_comm_(o.basis_comm_), linear_(o.linear_), basis_(o.basis_){
			static_assert(std::is_constructible<element_type, T>{}, "!");
		}

		auto skeleton() const {
			return skeleton_wrapper<field<basis_type, element_type>>(*this);
		}

		field(const field & coeff) = delete; // TODO make fields copyable
		field(field && coeff) = default;
		field & operator=(const field & coeff) = default;
		field & operator=(field && coeff) = default;

		//set to a scalar value
		field& operator=(element_type const& value){ // this makes sense only for zero?
			linear_.fill(value);
			return *this;
		}

		template<typename OtherT>
		field& operator=(field<basis_type, OtherT> const& o){
			static_assert( std::is_assignable<element_type&, OtherT>{}, "!" );
			assert( o.basis_ == basis_ and o.basis_comm_ == basis_comm_ );
			linear() = o.linear();
			return *this;
		}

		auto size() const {
			return basis_.size();
		}
		
		const auto & basis() const {
			return basis_;
		}

		auto cubic() const {
			return linear_.partitioned(basis_.cubic_dist(1).local_size()*basis_.cubic_dist(0).local_size()).partitioned(basis_.cubic_dist(0).local_size());
		}

		auto cubic() {
			return linear_.partitioned(basis_.cubic_dist(1).local_size()*basis_.cubic_dist(0).local_size()).partitioned(basis_.cubic_dist(0).local_size());
		}

		auto data() {
			return linear_.data();
		}

		auto data() const {
			return linear_.data();
		}

		auto & linear() const {
			return linear_;
		}

		auto & linear() {
			return linear_;
		}

		template <int dir = 2>
		friend void print_debug(const field & fld, const std::string & filename){

			std::ofstream file(filename);

			std::array<int, 3> point = {0, 0, 0};

			auto size = std::get<dir>(sizes(fld.cubic()));
			
			for(int ii = 0; ii < size; ii++){
				auto ip = ii + size/2;
				if(ip >= size) ip -= size;
				point[dir] = ip;
				auto rr = fld.basis().rvector(point);
				tfm::format(file, "%f %e %e\n", rr[dir], ::real(fld.cubic()[point[0]][point[1]][point[2]]), imag(fld.cubic()[point[0]][point[1]][point[2]]));
			}
		}

		auto & basis_comm() const {
			return basis_comm_;
		}


		auto complex() const {
			return field<basis::real_space, std::complex<element_type>>(*this);
		}

		field<basis::real_space, double> real() const {
			field<basis::real_space, double> real_field(skeleton());

			// Multi should be able to do this, but it causes a lot of compilation troubles
			//			
			real_field.linear() = boost::multi::blas::real(linear());
			
			//DATAOPERATIONS GPU::RUN 1D
		//	gpu::run(basis().part().local_size(),
		//					 [rp = begin(real_field.linear()), cp = begin(linear())] GPU_LAMBDA (auto ii){
		//						 rp[ii] = ::real(cp[ii]);
		//					 });
			return real_field;
		}
		
	private:
		mutable boost::mpi3::communicator basis_comm_;
		internal_array_type linear_;
		basis_type basis_;

	};
	
}

#if (not __INCLUDE_LEVEL__) or defined(UNIT_TEST) or defined(_TEST_BASIS_FIELD)
#if (not __INCLUDE_LEVEL__)
#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
int main( int argc, char* argv[] ) {
	boost::mpi3::environment env(argc, argv);
	return Catch::Session().run( argc, argv );
}
#endif

#include <basis/real_space.hpp>

#include <ions/unitcell.hpp>


TEST_CASE("Class basis::field", "[basis::field]"){

	using namespace Catch::literals;
	using math::vec3d;

	double ecut = 40.0;

	auto comm = boost::mpi3::environment::get_world_instance();
	
	ions::UnitCell cell(vec3d(10.0, 0.0, 0.0), vec3d(0.0, 4.0, 0.0), vec3d(0.0, 0.0, 7.0));
	basis::real_space rs(cell, input::basis::cutoff_energy(ecut), comm);

	basis::field<basis::real_space, double> ff(rs, comm);

	basis::field<basis::real_space, std::complex<double> > ff_complex = ff.complex();
	basis::field<basis::real_space, double> ff2 = ff_complex.real();

	ff2 = 0.;

	CHECK(( sizes(rs) == decltype(sizes(rs)){28, 11, 20} ));

	if(comm.size() == 1) CHECK(std::get<0>(sizes(ff.linear())) == 6160);
	if(comm.size() == 2) CHECK(std::get<0>(sizes(ff.linear())) == 3080);
	if(comm.size() == 4) CHECK(std::get<0>(sizes(ff.linear())) == 1540);

	if(comm.size() == 1) CHECK(std::get<0>(sizes(ff.cubic())) == 28);
	if(comm.size() == 2) CHECK(std::get<0>(sizes(ff.cubic())) == 14);
	if(comm.size() == 4) CHECK(std::get<0>(sizes(ff.cubic())) == 7);
	CHECK(std::get<1>(sizes(ff.cubic())) == 11);
	CHECK(std::get<2>(sizes(ff.cubic())) == 20);

	ff = 12.2244;

	for(int ii = 0; ii < rs.part().local_size(); ii++) CHECK(ff.linear()[ii] == 12.2244_a);	

	basis::field<basis::real_space, double> ff_copy(ff.skeleton());

	CHECK(std::get<1>(sizes(ff_copy.cubic())) == 11);
	CHECK(std::get<2>(sizes(ff_copy.cubic())) == 20);

	auto zff = ff.complex();
	
	static_assert(std::is_same<decltype(zff), basis::field<basis::real_space, std::complex<double>>>::value, "complex() should return a complex field");
	
	CHECK(std::get<1>(sizes(zff.cubic())) == 11);
	CHECK(std::get<2>(sizes(zff.cubic())) == 20);

	auto dff = zff.real();

	static_assert(std::is_same<decltype(dff), basis::field<basis::real_space, double>>::value, "real() should return a double field");

	CHECK(std::get<1>(sizes(dff.cubic())) == 11);
	CHECK(std::get<2>(sizes(dff.cubic())) == 20);

}

#endif

#endif
