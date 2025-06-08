cmake_minimum_required(VERSION 3.10)

project(Incomp2D_Project LANGUAGES Fortran)

# Set default build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# Fortran compiler configuration
set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O")

set(SOURCES
    driver.f90
    axycalc.f90
    hxycalc.f90
    output.f90
    press.f90
    qhat.f90
    sy.f90
    trigon.f90
    ufinal.f90
    fftcos.f90
    fftpackt.f90
)

add_executable(incomp2d_solver ${SOURCES})


