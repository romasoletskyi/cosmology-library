set(LIB_DIR "${cosmology_SOURCE_DIR}/include/cosmology")
set(HEADER_LIST "${LIB_DIR}/odesolver.h"
                "${LIB_DIR}/parameters.h"
                "${LIB_DIR}/perturbations.h"
                "${LIB_DIR}/utility.h")

set(SRC_LIST "perturbations.cpp")

add_library(cosmology_lib ${HEADER_LIST} ${SRC_LIST})
target_include_directories(cosmology_lib PUBLIC ../include)
target_link_libraries(cosmology_lib Eigen3::Eigen)
set_target_properties(cosmology_lib PROPERTIES LINKER_LANGUAGE CXX)
