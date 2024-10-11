
set(HEADER_LIST "${PROJECT_SOURCE_DIR}/include/ifd_lib/mesh.hpp")
add_library(ifd_lib mesh.cpp ${HEADER_LIST})

target_include_directories(ifd_lib PUBLIC ../include)

source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
