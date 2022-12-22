cmake_minimum_required(VERSION 3.23)
project(cosmology
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 3.3 REQUIRED)
find_package(GTest REQUIRED)

add_subdirectory(src)
add_subdirectory(apps)
add_subdirectory(tests)
