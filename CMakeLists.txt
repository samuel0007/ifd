cmake_minimum_required(VERSION 3.19)

project(
        IFD
        VERSION 1.0
        LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)


# include(FetchContent)
# FetchContent_Declare(matplotplusplus
#         GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
#         GIT_TAG origin/master)

#         FetchContent_GetProperties(matplotplusplus)
# if(NOT matplotplusplus_POPULATED)
#     FetchContent_Populate(matplotplusplus)
#     add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
# endif()

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)


add_subdirectory(src)

include_directories( ${Boost_INCLUDE_DIR} )


add_executable(main main.cpp)
# target_link_libraries(main PRIVATE ifd_lib Eigen3::Eigen ${Boost_LIBRARIES} matplot)
target_link_libraries(main PRIVATE ifd_lib Eigen3::Eigen ${Boost_LIBRARIES})

