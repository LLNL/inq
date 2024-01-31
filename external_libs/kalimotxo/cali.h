/* -*- indent-tabs-mode: t -*- */

#ifndef KALIMOTXO__CALI
#define KALIMOTXO__CALI

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <chrono>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <list>
#include <cassert>
#include <vector>
#include <algorithm>

#ifdef ENABLE_CALIPER

#include <caliper/cali.h>
#include <caliper/cali-manager.h>

#else

#define CALI_MARK_BEGIN(xx) kalimotxo::timer::mark_begin(xx);
#define CALI_MARK_END(xx) kalimotxo::timer::mark_end(xx);
#define CALI_CXX_MARK_SCOPE(xx) kalimotxo::timer __kalimotxo_scope_timer##__LINE__(xx);
#define CALI_CXX_MARK_FUNCTION kalimotxo::timer __kalimotxo_function_timer##__LINE__(__func__);
#define CALI_CXX_MARK_LOOP_BEGIN(xx, yy)
#define CALI_CXX_MARK_LOOP_END(xx)


namespace kalimotxo {

struct timer;

struct profile {
	profile():
		exc_time(0.0),
		inc_time(0.0),
		active_timer(NULL),
		call_count(0)
	{		
	}

	double exc_time;
	double inc_time;
	timer * active_timer;
	int call_count;
	
};

struct global_object {

	static auto & get() {
		static global_object go;
		return go;
	}
	
	global_object():
		enabled(false){
	}

	bool enabled;
	std::unordered_map<std::string, profile> timers_;
	std::list<decltype(timers_.begin())> stack_;
	
};


struct timer {

	std::string label_;
	std::chrono::time_point<std::chrono::high_resolution_clock> entry_;

	timer(std::string const & label):
		label_(label),
		entry_(std::chrono::high_resolution_clock::now())
	{
		if(not global_object::get().enabled) return;
		
		auto it = global_object::get().timers_.try_emplace(label_);
		global_object::get().stack_.push_back(it.first);
	}

	~timer(){
		if(not global_object::get().enabled) return;
		
		std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now() - entry_;
		global_object::get().timers_[label_].inc_time += elapsed_seconds.count();
		global_object::get().timers_[label_].exc_time += elapsed_seconds.count();
		global_object::get().timers_[label_].call_count++;
		global_object::get().stack_.pop_back();

		if(not global_object::get().stack_.empty()){
			global_object::get().stack_.back()->second.exc_time -= elapsed_seconds.count();
		}
	}

	static void mark_begin(std::string const & label){		

		if(not global_object::get().enabled) return;
		
		timer * tim = new timer(label);
		global_object::get().stack_.back()->second.active_timer = tim;
		
	}
	
	static void mark_end(std::string const & label){
		if(not global_object::get().enabled) return;

		auto active_prof = global_object::get().stack_.back();

		assert(active_prof->first == label);
		assert(active_prof->second.active_timer != NULL);

		delete active_prof->second.active_timer;
		active_prof->second.active_timer = NULL;
		
	}
	
};

}

namespace cali {

class ConfigManager {
	
public:
	
	void add(std::string const &) {
	}
	
	void start() {
		kalimotxo::global_object::get().enabled = true;
	}
	
	void flush() {

		assert(kalimotxo::global_object::get().stack_.empty());

		auto & timers = kalimotxo::global_object::get().timers_;

		std::vector<decltype(timers.begin())> sorted_timers;

		for(auto timer = timers.begin(); timer != timers.end(); ++timer) sorted_timers.push_back(timer);

		std::sort(sorted_timers.begin(), sorted_timers.end(), [](auto aa, auto bb) { return aa->second.exc_time > bb->second.exc_time; });

		auto file = fopen("profile.dat", "w");
		for(auto timer = sorted_timers.begin(); timer != sorted_timers.end(); ++timer){		
			fprintf(file, "  %-50s\t%15d\t%10.2f\t%10.2f\n", (*timer)->first.substr(0, 50).c_str(), (*timer)->second.call_count, (*timer)->second.inc_time, (*timer)->second.exc_time);
		}
		fclose(file);
	}
  
};

}

#endif

#endif
