cmake_minimum_required(VERSION 3.19)

project(
    IFD
    VERSION 1.0
    LANGUAGES CXX
)

find_package(Eigen3 REQUIRED)

add_subdirectory(src)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE ifd_lib)