cmake_minimum_required(VERSION 4.1)
get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(empty_spheres_reconstruction)
#set(CMAKE_BUILD_TYPE debug)
set(CMAKE_BUILD_TYPE release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
set(EIGEN_BUILD_DOC OFF)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
set( OFF)
FetchContent_MakeAvailable(Eigen)

find_package(CGAL REQUIRED COMPONENTS Core)

add_executable(${PROJECT_NAME} empty_spheres_reconstruction_3.cpp)
target_link_libraries(${PROJECT_NAME} CGAL::CGAL Eigen3::Eigen)
