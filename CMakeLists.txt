cmake_minimum_required (VERSION 2.6)

#set(CMAKE_C_COMPILER cl)
#set(CMAKE_CXX_COMPILER cl)

#set(CMAKE_INCLUDE_DIRECTORIES "E:\'Program Files (x86)'\IntelSWTools\compilers_and_libraries_2018.1.156\windows\compiler\include")

if(UNIX)
	set(CMAKE_Fortran_COMPILER ifort)
endif(UNIX)

project(Diploma)

enable_language(Fortran)

add_executable(Ogorod main.f90 Lsq.f90 healpix.f90 Galactic.f90 DataTypes.f90 ogorod.f90)

#find_package(ifqwin)

add_executable(Graphics Prog.f90 Lsq.f90 healpix.f90 Galactic.f90 DataTypes.f90 ogorod.f90 PROJECTION.F90 ColorPalette.f90)

set_target_properties(Graphics PROPERTIES COMPILE_OPTIONS "/libs:qwins" LINK_FLAGS "/subsystem:windows")

target_link_libraries(Graphics ifqwin)
