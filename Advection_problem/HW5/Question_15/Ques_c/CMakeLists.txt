# -----------------------------------------------------------------------
# Author: Bojun Zhang
# Date  : 2025-04-15
# File  : CMakeLists.txt
# -----------------------------------------------------------------------

cmake_minimum_required(VERSION 3.10)

# Project name and language
project(AdvectionSolver LANGUAGES Fortran)

# Set Fortran standard
set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED True)

# List your source files
set(SOURCES
    Advection.f90
    find_slope.f90
)

# Define the executable
add_executable(advection_solver ${SOURCES})

add_custom_target(clean_run_files
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/output.dat
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/fort.*
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/*.txt
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/*.pdf
    COMMENT "Cleaning runtime-generated files"
)


# Optional: set compile flags
target_compile_options(advection_solver PRIVATE -O2 -Wall)

# Optional: set output directory
set_target_properties(advection_solver PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Message for build info
message(STATUS "Executable will be placed in ${CMAKE_BINARY_DIR}")

