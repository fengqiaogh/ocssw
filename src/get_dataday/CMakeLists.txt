cmake_minimum_required(VERSION 3.12)

#
# Make a version file containing the current version from git.
#
set (VERSION_MAJOR 1)
set (VERSION_MINOR 1)
set (VERSION_PATCH 0)
include(GetOcsswVersion)

# set(Boost_DEBUG 1)
find_package(Boost REQUIRED)

find_package(NetCDF REQUIRED COMPONENTS CXX)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})
list(APPEND LIBS timeutils)
list(APPEND LIBS piutils)


add_executable(get_dataday
  get_dataday.cpp
)

target_link_libraries(get_dataday ${LIBS})

target_include_directories(get_dataday SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})

install (TARGETS get_dataday DESTINATION bin)

