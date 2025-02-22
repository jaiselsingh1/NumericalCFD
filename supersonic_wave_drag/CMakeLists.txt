cmake_minimum_required(VERSION 3.22)
project(hw)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -Wpedantic -Wextra")

# External Libraries
get_filename_component(PARENT_DIR ${PROJECT_SOURCE_DIR} DIRECTORY)
list(APPEND CMAKE_PREFIX_PATH "${PARENT_DIR}/third_party/matplotplusplus/install")
find_package(Matplot++ REQUIRED)
list(APPEND CMAKE_PREFIX_PATH "${PARENT_DIR}/third_party/eigen/install")
find_package(Eigen3 REQUIRED)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/source
)

file(GLOB all_SRCS
    "${PROJECT_SOURCE_DIR}/include/*.h"
    "${PROJECT_SOURCE_DIR}/include/*.hpp"
    "${PROJECT_SOURCE_DIR}/source/*.cpp"
    "${PROJECT_SOURCE_DIR}/source/*.c"
)

add_executable(${CMAKE_PROJECT_NAME} ${all_SRCS})
target_link_libraries(${CMAKE_PROJECT_NAME} m)
target_link_libraries(${CMAKE_PROJECT_NAME} Matplot++::matplot)
target_link_libraries(${CMAKE_PROJECT_NAME} Eigen3::Eigen)