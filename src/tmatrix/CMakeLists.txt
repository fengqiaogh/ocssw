cmake_minimum_required(VERSION 3.12)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(LIBS genutils)

find_package(NetCDF REQUIRED COMPONENTS C CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

add_executable(tmatrix
    tm_tmd.lp.f
    lpd.f    
    tmatrix.cpp
    tmatrix_input.cpp
    TmParamsReader.cpp
    TmProcess.cpp
)

target_link_libraries(tmatrix ${LIBS})

install (TARGETS tmatrix DESTINATION bin)
