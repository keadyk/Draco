#-----------------------------*-cmake-*----------------------------------------#
# file   compton/CMakeLists.txt
# author Kendra Keady
# date   2016 Oct 14
# brief  Generate build project files for compton.
# note   Copyright (C) 2010-2016, Los Alamos National Security, LLC.
#        All rights reserved.
#------------------------------------------------------------------------------#
# $Id$
#------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.0.0)
project( compton CXX )

# ---------------------------------------------------------------------------- #
# Source files
# ---------------------------------------------------------------------------- #

file( GLOB template_implementations *.t.hh *.i.hh )
file( GLOB sources *.cc )
file( GLOB explicit_instantiations *_pt.cc )
file( GLOB headers *.hh )

# Make the header files available in the IDE.
if( MSVC_IDE OR ${CMAKE_GENERATOR} MATCHES Xcode )
   list( APPEND sources ${headers} )
endif()

# ---------------------------------------------------------------------------- #
# Directories to search for include directives
# ---------------------------------------------------------------------------- #

include_directories(
   ${PROJECT_SOURCE_DIR}       # sources
   ${draco_src_dir_SOURCE_DIR} # ds++ header files
   ${dsxx_BINARY_DIR}          # ds++/config.h
   )

# ---------------------------------------------------------------------------- #
# Build package library
# ---------------------------------------------------------------------------- #

add_component_library(
   TARGET       Lib_compton
   TARGET_DEPS  "Lib_units;GSL::gsl"
   LIBRARY_NAME ${PROJECT_NAME}
   SOURCES      "${sources}"
)

# ---------------------------------------------------------------------------- #
# Installation instructions
# ---------------------------------------------------------------------------- #

install( TARGETS Lib_compton EXPORT draco-targets DESTINATION ${DBSCFGDIR}lib )
install( FILES ${headers} ${template_implementations} 
         DESTINATION ${DBSCFGDIR}include/compton )

# ---------------------------------------------------------------------------- #
# Unit tests
# ---------------------------------------------------------------------------- #

if( BUILD_TESTING )
   add_subdirectory( test )
endif()

# ---------------------------------------------------------------------------- #
# Autodoc
# ---------------------------------------------------------------------------- #

process_autodoc_pages()
