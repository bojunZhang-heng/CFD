cmake_minimum_required(VERSION 3.10)
project(incomp2d Fortran)

set(CMAKE_Fortran_COMPILER gfortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -Wall")

# 所有源文件（自动列出当前目录下所有 .f90 文件）
file(GLOB SRC_FILES *.f90)

# 设置模块输出目录
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# 构建可执行文件
add_executable(incomp2d ${SRC_FILES})
set_target_properties(incomp2d PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

