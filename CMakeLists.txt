cmake_minimum_required(VERSION 3.23)
project(cosmology
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_ASAN "-g -fsanitize=address,undefined -fno-sanitize-recover=all"
        CACHE STRING "Compiler flags in asan build"
        FORCE)

include(FetchContent)
FetchContent_Declare(complex_bessel
        GIT_REPOSITORY "https://github.com/joeydumont/complex_bessel.git"
        GIT_TAG 1c0b9ea839077031f4b53e77a23a4f5e3e6f84fb
        )
FetchContent_MakeAvailable(complex_bessel)

find_package(Eigen3 3.3 REQUIRED)
find_package(GTest REQUIRED)

add_subdirectory(src)
add_subdirectory(apps)
add_subdirectory(tests)
