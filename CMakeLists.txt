cmake_minimum_required(VERSION 3.15)
project(cma_laboratory_work_2)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_executable(cma_laboratory_work_2 main.cpp tasks/qr_algorithm.h tasks/constants.h tasks/constants.cpp tasks/matrix.h tasks/matrix.cpp)