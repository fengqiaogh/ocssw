cmake_minimum_required(VERSION 3.12)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(LIBS l1 seawifs dfutils)

#find_package(NetCDF REQUIRED)
#include_directories(${NETCDF_INCLUDE_DIRS})
#list(APPEND LIBS ${NETCDF_LIBRARIES})

ADD_LIBRARY(l2
  readL2scan.c
  setupflags.c
  get_product_table.c
  expand3D.cpp
)

target_link_libraries(l2 ${LIBS})
target_include_directories(l2 PUBLIC ".")

install (TARGETS l2 DESTINATION lib)
