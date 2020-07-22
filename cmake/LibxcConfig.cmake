# LibxcConfig.cmake
# ------------------
#
# Libxc cmake module.
# This module sets the following variables in your project:
#
# ::
#
#   Libxc_FOUND - true if Libxc and all required components found on the system
#   Libxc_VERSION - Libxc version in format Major.Minor.Release
#   Libxc_INCLUDE_DIRS - Directory where Libxc header is located.
#   Libxc_INCLUDE_DIR - same as DIRS
#   Libxc_DEFINITIONS: Definitions necessary to use Libxc, namely USING_Libxc.
#   Libxc_LIBRARIES - Libxc library to link against.
#   Libxc_LIBRARY - same as LIBRARIES
#
#
# Available components: shared static
#
# ::
#
#   shared - search for only shared library
#   static - search for only static library
#
#
# Exported targets:
#
# ::
#
# If Libxc is found, this module defines the following :prop_tgt:`IMPORTED`
# target. Target is shared _or_ static, so, for both, use separate, not
# overlapping, installations. ::
#
#   Libxc::xc - the main Libxc library with header & defs attached.
#
#
# Suggested usage:
#
# ::
#
#   find_package(Libxc)
#   find_package(Libxc 3.0.0 EXACT CONFIG REQUIRED COMPONENTS shared)
#
#
# The following variables can be set to guide the search for this package:
#
# ::
#
#   Libxc_DIR - CMake variable, set to directory containing this Config file
#   CMAKE_PREFIX_PATH - CMake variable, set to root directory of this package
#   PATH - environment variable, set to bin directory of this package
#   CMAKE_DISABLE_FIND_PACKAGE_Libxc - CMake variable, disables
#     find_package(Libxc) when not REQUIRED, perhaps to force internal build


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was LibxcConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

set(PN Libxc)
set (_valid_components
    static
    shared
)

# find includes
unset(_temp_h CACHE)
find_path(_temp_h
          NAMES //xc.h
          PATHS ${PACKAGE_PREFIX_DIR}/include
          NO_DEFAULT_PATH)
if(_temp_h)
    set(${PN}_INCLUDE_DIR "${_temp_h}")
    set(${PN}_INCLUDE_DIRS ${${PN}_INCLUDE_DIR})
else()
    set(${PN}_FOUND 0)
    if(NOT CMAKE_REQUIRED_QUIET)
        message(STATUS "${PN}Config missing component: header (${PN}: ${_temp_h})")
    endif()
endif()

# find library: shared, static, or whichever
set(_hold_library_suffixes ${CMAKE_FIND_LIBRARY_SUFFIXES})
list(FIND ${PN}_FIND_COMPONENTS "shared" _seek_shared)
list(FIND ${PN}_FIND_COMPONENTS "static" _seek_static)
if(_seek_shared GREATER -1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX})
elseif(_seek_static GREATER -1)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()
unset(_temp CACHE)
find_library(_temp
             NAMES xc
             PATHS ${PACKAGE_PREFIX_DIR}/lib
             NO_DEFAULT_PATH)
if(_temp)
    set(${PN}_LIBRARY "${_temp}")
    if(_seek_shared GREATER -1)
        set(${PN}_shared_FOUND 1)
    elseif(_seek_static GREATER -1)
        set(${PN}_static_FOUND 1)
    endif()
else()
    if(_seek_shared GREATER -1)
        if(NOT CMAKE_REQUIRED_QUIET)
            message(STATUS "${PN}Config missing component: shared library (${PN}: ${_temp})")
        endif()
    elseif(_seek_static GREATER -1)
        if(NOT CMAKE_REQUIRED_QUIET)
            message(STATUS "${PN}Config missing component: static library (${PN}: ${_temp})")
        endif()
    else()
        set(${PN}_FOUND 0)
        if(NOT CMAKE_REQUIRED_QUIET)
            message(STATUS "${PN}Config missing component: library (${PN}: ${_temp})")
        endif()
    endif()
endif()
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_hold_library_suffixes})
set(${PN}_LIBRARIES ${${PN}_LIBRARY})
set(${PN}_DEFINITIONS USING_${PN})

check_required_components(${PN})

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already built this as a subproject
#-----------------------------------------------------------------------------
if(NOT TARGET ${PN}::xc)
    include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets.cmake")
endif()
