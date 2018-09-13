cmake_minimum_required(VERSION 3.12)
project(hw1_hpc)
FIND_PACKAGE(MPI REQUIRED)
set(CMAKE_CXX_STANDARD 14)


set(INCLUDE_FILES include/partition.h
                  include/quicksort.h )


set(SOURCE_FILES src/main.cpp
                 src/partition.cpp
                 src/quicksort.cpp)


add_executable(hw1_hpc ${INCLUDE_FILES} ${SOURCE_FILES})
target_include_directories(hw1_hpc PUBLIC include)