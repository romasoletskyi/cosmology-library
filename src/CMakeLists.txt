set(LIB_DIR "${cosmology_SOURCE_DIR}/include/cosmology")
set(HEADER_LIST "${LIB_DIR}/odesolver/odesolver.h"
                "${LIB_DIR}/odesolver/odesystem.h"
                "${LIB_DIR}/odesolver/parser.h"
                "${LIB_DIR}/odesolver/odewalker.h"
                "${LIB_DIR}/utility/math.h"
                "${LIB_DIR}/utility/graph.h"
                "${LIB_DIR}/utility/spline.h"
                "${LIB_DIR}/utility/utility.h"
                "${LIB_DIR}/parameters.h"
                "${LIB_DIR}/perturbations.h"
                "${LIB_DIR}/monte_carlo.h")

set(SRC_LIST "perturbations.cpp"
             "utility/spline.cpp")

add_library(cosmology_lib ${HEADER_LIST} ${SRC_LIST})
target_include_directories(cosmology_lib PUBLIC ../include)
target_link_libraries(cosmology_lib Eigen3::Eigen)
set_target_properties(cosmology_lib PROPERTIES LINKER_LANGUAGE CXX)
