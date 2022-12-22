add_executable(odesolver odesolver.cpp)
target_link_libraries(odesolver PRIVATE cosmology_lib GTest::gtest_main)