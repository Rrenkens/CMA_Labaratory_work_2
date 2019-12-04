cmake_minimum_required(VERSION 3.15)
project(cma_laboratory_work_2)

set(CMAKE_CXX_STANDARD 14)

add_executable(cma_laboratory_work_2 main.cpp tasks/qr_algorithm.h tasks/constants.h tasks/constants.cpp tasks/matrix.h tasks/matrix.cpp)