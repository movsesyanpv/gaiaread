cmake_minimum_required (VERSION 2.6)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /libs:qwin")

#set(CMAKE_INCLUDE_DIRECTORIES "E:\'Program Files (x86)'\IntelSWTools\compilers_and_libraries_2018.1.156\windows\compiler\include")

project(Diploma)

enable_language(Fortran)

add_executable(Ogorod main.f90 Lsq.f90 healpix.f90 Galactic.f90 DataTypes.f90 ogorod.f90)

#find_package(ifqwin)

add_executable(Graphics Prog.f90 Lsq.f90 healpix.f90 Galactic.f90 DataTypes.f90 ogorod.f90 PROJECTION.F90 ColorPalette.f90)

set_target_properties(Graphics PROPERTIES LINK_FLAGS "/subsystem:windows")

target_link_libraries(Graphics ifqwin)
