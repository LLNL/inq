/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INTERFACE__HISTORY_FILE
#define INQ__INTERFACE__HISTORY_FILE

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

namespace inq {
namespace interface {

struct {

	std::string entries_;
	
public:
	
	void add_entry(int argc, char ** argv) {

		for(int iarg = 0; iarg < argc; iarg++) {
			auto arg = std::string(argv[iarg]);
			auto & npos = std::string::npos;
			if(arg.find(' ') != npos or arg.find('(') != npos or arg.find(')') != npos){
				arg = '\"' + arg + '\"';
			}
			if(iarg > 0) entries_ += ' ';
			entries_ += arg;
		}
		entries_ += '\n';
	}

	void write() const {
		auto history_file = std::ofstream(".inq_history", std::ofstream::app);
		history_file << entries_;
	}
	
} history_file;

}
}

#endif

#ifdef INQ_INTERFACE_HISTORY_FILE_UNIT_TEST
#undef INQ_INTERFACE_HISTORY_FILE_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;

}
#endif
