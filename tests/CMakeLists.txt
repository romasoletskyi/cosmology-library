add_executable(odesolver odesolver.cpp)
add_executable(spline spline.cpp)
add_executable(math math.cpp)

target_link_libraries(odesolver PRIVATE cosmology_lib GTest::gtest_main)
target_link_libraries(spline PRIVATE cosmology_lib GTest::gtest_main)
target_link_libraries(math PRIVATE cosmology_lib GTest::gtest_main)