cmake_minimum_required(VERSION 3.23)
project(Conventional)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

set(INCLUDE include/perturbations.h
            include/odesolver.h
            include/utility.h
            include/parameters.h)

add_executable(odesolver ${INCLUDE} tests/odesolver.cpp)
target_link_libraries(odesolver ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES})

add_executable(perturbations ${INCLUDE} tests/perturbations.cpp)
