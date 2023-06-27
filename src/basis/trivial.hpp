/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__BASIS__TRIVIAL
#define INQ__BASIS__TRIVIAL

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cassert>
#include <array>

#include <basis/base.hpp>

namespace inq {
namespace basis {

	/*
		This is a class that implements a very simple basis object. Useful for testing.
	*/
	
  class trivial : public base {

  public:

		using reciprocal_space = void;
		
		const static int dimension = 1;
		
		trivial(const long size, parallel::communicator comm):
			base(size, comm),
      size_(size){
		}
				
		trivial(const long size, boost::mpi3::communicator comm):
			trivial(size, parallel::communicator{std::move(comm)}){
		}
		
    long size() const {
      return size_;
    }

		friend auto sizes(const trivial & ss){
      return std::array<long, 1>{ss.size_};
		}

		auto sizes() const {
			return std::array<long, 1>{size_};
		}

    double volume_element() const {
      return 1.0/size_;
    }    

    friend bool operator==(trivial const & b1, trivial const & b2){
      return b1.size_ == b2.size_;
    }
    
	protected:

    long size_;
		
  };

}
}
#endif

#ifdef INQ_BASIS_TRIVIAL_UNIT_TEST
#undef INQ_BASIS_TRIVIAL_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
  
	using namespace inq;
	using namespace Catch::literals;
  
}
#endif
