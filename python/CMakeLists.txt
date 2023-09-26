# -*- indent-tabs-mode: t -*-
# Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

########################
# PYBIND11
########################

FetchContent_Declare(
  pybind11
	GIT_REPOSITORY https://github.com/pybind/pybind11.git
	GIT_TAG        v2.11.1
  GIT_PROGRESS   TRUE
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(pybind11)

file (GLOB PYTHON_FILES "*.cpp" "*.hpp")

if(ENABLE_CUDA)
	set_source_files_properties(${PYTHON_FILES} PROPERTIES LANGUAGE CUDA)   
endif()

pybind11_add_module(pinq MODULE
	${PYTHON_FILES}
)

target_link_libraries(pinq PUBLIC)
target_link_libraries(pinq PRIVATE inq)

#DO NOT INSTALL FROM HERE
#
#if (NOT CMAKE_INSTALL_PYTHON_LIBDIR)
#  set(CMAKE_INSTALL_PYTHON_LIBDIR
#    "${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages")
#endif ()
#install(TARGETS pinq
#  COMPONENT python
#  LIBRARY DESTINATION "${CMAKE_INSTALL_PYTHON_LIBDIR}"
#)

