cmake_minimum_required(VERSION 3.23)
project(cosmology
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_ASAN "-g -fsanitize=address,undefined -fno-sanitize-recover=all"
        CACHE STRING "Compiler flags in asan build"
        FORCE)

find_package(Eigen3 3.3 REQUIRED)
find_package(GTest REQUIRED)

add_subdirectory(src)
add_subdirectory(apps)
add_subdirectory(tests)
