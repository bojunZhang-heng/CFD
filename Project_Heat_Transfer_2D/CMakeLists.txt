cmake_minimum_required(VERSION 3.10)
project(incomp2d Fortran)

SET(CMAKE_Fortran_COMPILER gfortran)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -Wall")

SET(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
SET(MODEL_DIR ${CMAKE_SOURCE_DIR}/Model_data)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# ==================================
# Include source 
# ==================================
SET(SRC_FILES
  ${SRC_DIR}/Driver.f90
  ${SRC_DIR}/fftcos.f90
  ${SRC_DIR}/fftpackt.f90
  ${SRC_DIR}/hxycalc.f90
  ${SRC_DIR}/axycalc.f90
  ${SRC_DIR}/qhat.f90
  ${SRC_DIR}/sy.f90
  ${SRC_DIR}/trigon.f90
  ${SRC_DIR}/ufinal.f90
  ${SRC_DIR}/press.f90
  ${SRC_DIR}/output.f90
  ${MODEL_DIR}/mod_fft_data.f90
  ${MODEL_DIR}/mod_param.f90
  ${MODEL_DIR}/mod_variables.f90

   )

add_executable(incomp2d ${SRC_FILES})

set_target_properties(incomp2d PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})
