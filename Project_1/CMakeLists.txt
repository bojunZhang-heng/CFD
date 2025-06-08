cmake_minimum_required(VERSION 3.10)
project(Project_1 Fortran)

# Set default build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# Fortran compiler configuration
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra")

# Different flags for Debug and Release
set(CMAKE_Fortran_FLAGS_DEBUG "-g -fcheck=all -fbacktrace")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")

# Source files
set(SOURCES
    Advection_1D.f90
    find_slope.f90
)

# Executable target
add_executable(Advection1D ${SOURCES})

# Installation (optional)
install(TARGETS Advection1D DESTINATION bin)
