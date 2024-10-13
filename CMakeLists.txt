cmake_minimum_required(VERSION 3.19)

project(
        IFD
        VERSION 1.0
        LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)


find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
add_subdirectory(src)

include_directories( ${Boost_INCLUDE_DIR} )

add_executable(main main.cpp)
target_link_libraries(main PRIVATE ifd_lib Eigen3::Eigen ${Boost_LIBRARIES})