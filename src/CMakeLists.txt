
set(HEADER_LIST "${PROJECT_SOURCE_DIR}/include/ifd_lib/mesh.hpp")
add_library(ifd_lib mesh.cpp ${HEADER_LIST})

target_include_directories(ifd_lib PUBLIC ../include)
include_directories(
        ${PYTHON3_INCLUDE_DIRS}
        "/usr/include/python3.8"
)
# target_link_libraries(ifd_lib PRIVATE Eigen3::Eigen matplot)
target_link_libraries(ifd_lib PRIVATE Eigen3::Eigen)


source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
