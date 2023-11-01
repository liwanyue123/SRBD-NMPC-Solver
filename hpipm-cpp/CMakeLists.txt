
#####################
## Build hpipm-cpp ##
#####################
# find blasfeo
find_package(blasfeo REQUIRED)
# find hpipm
find_package(hpipm REQUIRED)
# find Eigen
find_package(Eigen3 REQUIRED)

add_library(
  hpipm_cpp
  src/detail/d_ocp_qp_dim_wrapper.cpp
  src/detail/d_ocp_qp_wrapper.cpp
  src/detail/d_ocp_qp_sol_wrapper.cpp
  src/detail/d_ocp_qp_ipm_arg_wrapper.cpp
  src/detail/d_ocp_qp_ipm_ws_wrapper.cpp
  src/ocp_qp_dim.cpp
  src/ocp_qp_ipm_solver_settings.cpp
  src/ocp_qp_ipm_solver_statistics.cpp
  src/ocp_qp_ipm_solver.cpp
)
target_link_libraries(
  hpipm_cpp
  PUBLIC
  blasfeo
  hpipm
)
target_include_directories(
  hpipm_cpp
  PUBLIC
  ${EIGEN3_INCLUDE_DIR}
  $<BUILD_INTERFACE:${blasfeo_DIR}/../include>
  $<BUILD_INTERFACE:${hpipm_DIR}/../include>
  include
)




