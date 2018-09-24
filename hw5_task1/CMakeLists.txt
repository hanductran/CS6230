cmake_minimum_required(VERSION 3.12)
project(hw5_hpc)

set(CMAKE_CXX_STANDARD 14)
FIND_PACKAGE(MPI REQUIRED)
find_package(OpenMP REQUIRED)

if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif()

set(INCLUDE_FILES include/simple_bcast.h
                  include/log_bcast.h)

set(SOURCE_FILES src/main.cpp
                 src/simple_bcast.cpp
                 src/log_bcast.cpp)

add_executable(hw5_task1 ${INCLUDE_FILES} ${SOURCE_FILES})
target_include_directories(hw5_task1 PUBLIC include)
target_include_directories(hw5_task1 PRIVATE ${MPI_INCLUDE_PATH})
target_link_libraries(hw5_task1 ${MPI_LIBRARIES} m)