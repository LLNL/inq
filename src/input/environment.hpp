/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__INPUT__ENVIRONMENT
#define INQ__INPUT__ENVIRONMENT

/*
 Copyright (C) 2020 Xavier Andrade

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

#include <utils/merge_optional.hpp>

#include <mpi3/environment.hpp>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>

#include <optional>
#include <cassert>

namespace inq {
namespace input {

  class environment {

  public:
    
    environment(int argc, char** argv):
      mpi_env_(argc, argv)
    {
			if(mpi_env_.get_world_instance().rank() == 0){
				calimgr_.add("runtime-report");
				calimgr_.start();
			}

			CALI_MARK_BEGIN("inq_environment");
			
    }

		~environment(){

			CALI_MARK_END("inq_environment");
			
			if(mpi_env_.get_world_instance().rank() == 0){
				calimgr_.flush(); // write performance results
			}
		}
      

  private:
    boost::mpi3::environment mpi_env_;
		cali::ConfigManager calimgr_;
		
  };
    
}
}

////////////////////////////////////////////////////////

#ifdef INQ_INPUT_ENVIRONMENT_UNIT_TEST
#undef INQ_INPUT_ENVIRONMENT_UNIT_TEST

#endif
   
#endif
