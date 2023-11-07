cmake_minimum_required(VERSION 3.5.1)
project(stance_nmpc_test)
# set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

# add_compile_options(-O2 -g -mavx -DNDEBUG -fPIC -std=gnu++17)
# add_compile_options(-mfma -mavx -O2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_BUILD_TYPE Release)



add_subdirectory(hpipm-cpp)
add_subdirectory(dynamic)
include_directories(
  include
)

find_package(yaml-cpp REQUIRED)
# find_package(matplotlib_cpp)
# find_package(PythonLibs REQUIRED)
# include_directories(${PYTHON_INCLUDE_DIRS})


# find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)


add_executable(${PROJECT_NAME}
  NMPC_solver.cpp
)

target_link_libraries(${PROJECT_NAME} hpipm_cpp flow_dynamic yaml-cpp  )
