# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required (VERSION 3.20 FATAL_ERROR)

# Set the project name
project (HOMEDF)

# detect operating system and host processor
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")

# Select flags.
SET(CMAKE_CXX_FLAGS "-fopenmp -O3 -DNDEBUG -Wall -march=native")
#SET(CMAKE_CXX_FLAGS "-fopenmp -g")

# try conditional compilation
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)

# check results and add flag
if(COMPILER_SUPPORTS_CXX17)
    set(CMAKE_CXX_STANDARD 17)
    message(STATUS "Using the compiler C++17 ${CMAKE_CXX_COMPILER}")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++17 support. Please use a different C++ compiler.")
endif()

##################################################################################
#VTK
find_package(VTK REQUIRED)
include_directories(SYSTEM ${VTK_INCLUDE_DIRS})
link_libraries(${VTK_LIBRARIES})

##################################################################################
#Eigen
find_package(Eigen3 REQUIRED)
if(EIGEN3_INCLUDE_DIR)
  include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
endif()

##################################################################################
#MPI
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
link_libraries(${MPI_CXX_LIBRARIES})

##################################################################################
#gtest
find_package(GTest REQUIRED)
include_directories(SYSTEM ${GTEST_INCLUDE_DIR})
link_libraries(GTest::gtest)
#link_libraries(GTest::GTest)
#link_libraries(${GTEST_LIBRARIES})
##################################################################################
#adding executable for the main
add_executable (Main Main_files/Solve_F.cpp)
add_executable (Closing_Flux Main_files/S5.cpp)

#adding executable for the tests
add_executable (Test Test/TEST_main.cpp)


message(STATUS "Generated build with CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
