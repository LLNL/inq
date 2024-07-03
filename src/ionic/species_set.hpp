/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__IONIC__SPECIES_SET
#define INQ__IONIC__SPECIES_SET

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <ionic/species.hpp>

namespace inq {
namespace ionic {

class species_set {

	using container_type = 	std::unordered_map<std::string, ionic::species>;

	container_type list_;

public:
	
	auto size() const {
		return (long) list_.size();
	}

	auto & list() const {
		return list_;
	}

	auto insert(ionic::species const & sp) {
		list_.insert_or_assign(sp.symbol(), sp);
	}

	auto contains(std::string const & symbol) const{
		return list_.find(symbol) != list_.end();
	}

	auto & operator[](std::string const & symbol) const{
		return list_.at(symbol);
	}

	struct const_iterator {

		container_type::const_iterator base_iter;
		
		auto operator!=(const_iterator const & other) const {
			return base_iter != other.base_iter;
		}

		auto operator++() {
			return const_iterator{++base_iter};
		}

		auto operator*() const {
			return base_iter->second;
		}

	};
	
	auto begin() const {
		return const_iterator{list_.begin()};
	}

	auto end() const {
		return const_iterator{list_.end()};
	}
	
};

}
}
#endif

#ifdef INQ_IONIC_SPECIES_SET_UNIT_TEST
#undef INQ_IONIC_SPECIES_SET_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {
	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
	
	auto list = ionic::species_set();

	list.insert("C");

	CHECK(list.contains("C"));
	CHECK(not list.contains("F"));

	auto local_h = inq::ionic::species("H").symbol("Hloc").pseudo(inq::config::path::unit_tests_data() + "H.blyp-vbc.UPF");
	list.insert(local_h);

	CHECK(list.contains("Hloc"));
	CHECK(not list.contains("H"));
	CHECK(list.size() == 2);
	CHECK(list["Hloc"].mass() == 1837.17994584_a);
	CHECK(list["Hloc"].file_path() == inq::config::path::unit_tests_data() + "H.blyp-vbc.UPF");

	list.insert("U");

	CHECK(list.size() == 3);

	auto count = 0;
	auto has_h = false;
	auto has_c = false;
	auto has_u = false;
	for(auto const & species : list) {
		if(species.symbol() == "Hloc") has_h = true;
		if(species.symbol() == "C") has_c = true;
		if(species.symbol() == "U") has_u = true;
		count++;
	}
	CHECK(has_h);
	CHECK(has_c);
	CHECK(has_u);
	CHECK(count == 3);
	
}

#endif
