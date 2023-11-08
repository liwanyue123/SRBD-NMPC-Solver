cmake_minimum_required(VERSION 3.5.1)
project(SRBD_NMPC)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

add_compile_options(-O3 -mfma -mavx -fPIC -std=gnu++17)
set(CMAKE_BUILD_TYPE Release)

add_subdirectory(hpipm-cpp)
add_subdirectory(dynamics)
include_directories(
  include
)

find_package(yaml-cpp REQUIRED)
 
add_executable(${PROJECT_NAME}
  NMPC_solver.cpp
)

target_link_libraries(${PROJECT_NAME} hpipm_cpp SRBD_model yaml-cpp  )
