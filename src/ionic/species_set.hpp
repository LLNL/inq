/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__IONIC__SPECIES_SET
#define INQ__IONIC__SPECIES_SET

// Copyright (C) 2019-2024 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <ionic/species.hpp>
#include <pseudopod/set_id.hpp>

namespace inq {
namespace ionic {

class species_set {

	using container_type = 	std::unordered_map<std::string, ionic::species>;

	container_type list_;
	pseudo::set_id pseudo_set_;

public:

	species_set(pseudo::set_id const & pseudos = pseudo::set_id::pseudodojo_pbe()):
		pseudo_set_(pseudos)
	{
	}

	auto & pseudopotentials() const {
		return pseudo_set_;
	}

	auto & pseudopotentials() {
		return pseudo_set_;
	}

	auto size() const {
		return (long) list_.size();
	}

	auto & list() const {
		return list_;
	}

	auto insert(ionic::species const & sp) {
		list_.insert_or_assign(sp.symbol(), sp);
	}

	auto contains(std::string const & symbol) const {
		return list_.find(symbol) != list_.end();
	}

	auto & operator[](std::string const & symbol) const {
		return list_.at(symbol);
	}

	auto & operator[](std::string const & symbol) {
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

	void save(parallel::communicator & comm, std::string const & dirname) const {
		auto error_message = "INQ error: Cannot save species_set to directory '" + dirname + "'.";

		utils::create_directory(comm, dirname);
		
		utils::save_value    (comm, dirname + "/pseudo_set", pseudo_set_, error_message);
		utils::save_container(comm, dirname + "/symbols",    list_,       error_message, [](auto el) { return el.first; });

		for(auto const & species : *this) species.save(comm, dirname + "/species_" + species.symbol());
	}
	
	static auto load(std::string const & dirname) {
		auto error_message = "INQ error: Cannot load species_set from directory '" + dirname + "'.";

		auto read = ionic::species_set();
		
		utils::load_value(dirname + "/pseudo_set", read.pseudo_set_, error_message);

		auto symbols = std::vector<std::string>{};
		utils::load_vector(dirname + "/symbols", symbols);

		for(auto const & sym : symbols) {
			auto species = ionic::species::load(dirname + "/species_" + sym);
			read.insert(species);
		}

		return read;
	}

	template<class OStream>
	friend OStream & operator<<(OStream & out, species_set const & self){

		using namespace magnitude;

		out << "Species (" << self.size() << "):\n";

		for(auto const & species : self) {
			out << "  " << species;
			if(species.default_pseudo()) out << " *";
			out << std::endl;
		}

		out << "\n  * Uses default pseudopotential set: '" << self.pseudopotentials() << "'" << std::endl;
		
		return out;
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

	parallel::communicator comm{boost::mpi3::environment::get_world_instance()};
	
	auto list = ionic::species_set();

	list.pseudopotentials() = pseudo::set_id::sg15();

	CHECK(list.pseudopotentials() == pseudo::set_id::sg15());
	
	list.insert("C");

	CHECK(list.contains("C"));
	CHECK(not list.contains("F"));

	auto local_h = ionic::species("H").symbol("Hloc").pseudo_file(config::path::unit_tests_data() + "H.blyp-vbc.UPF");
	list.insert(local_h);

	CHECK(list.contains("Hloc"));
	CHECK(not list.contains("H"));
	CHECK(list.size() == 2);
	CHECK(list["Hloc"].mass() == 1837.17994584_a);
	CHECK(list["Hloc"].file_path() == config::path::unit_tests_data() + "H.blyp-vbc.UPF");

	list.insert(ionic::species("U").mass(235));

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

	list.save(comm, "save_species_set");
	auto read_list = list.load("save_species_set");

	CHECK(read_list.size() == 3);
	CHECK(read_list.pseudopotentials() == pseudo::set_id::sg15());
	CHECK(read_list.contains("C"));
	CHECK(read_list.contains("Hloc"));
	CHECK(not read_list.contains("H"));
	CHECK(read_list["Hloc"].mass() == 1837.17994584_a);
	CHECK(read_list["Hloc"].file_path() == config::path::unit_tests_data() + "H.blyp-vbc.UPF");
	CHECK(read_list.contains("U"));

	list.insert(ionic::species("Te").nofilter());
	list.insert(ionic::species("Vn").symbol("Vn_ccecp").pseudo_set(pseudo::set_id::ccecp()).nofilter());

	std::cout << list << std::endl;
	
}

#endif
