cmake_minimum_required(VERSION 3.12)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(LIBS l2 nav)

set (NETCDF_F90 "YES")
find_package(NetCDF REQUIRED COMPONENTS C CXX Fortran)
include_directories(${NETCDF_INCLUDE_DIRS})
list(APPEND LIBS ${NETCDF_LIBRARIES})

find_package(Levmar REQUIRED)
include_directories(${Levmar_INCLUDE_DIRS})
list(APPEND LIBS ${Levmar_LIBRARIES})
  
find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIRS})
list(APPEND LIBS ${LAPACK_LIBRARIES})

add_executable(dtdb
    dtdb.cpp
    DDProcess.cpp
    DDAlgorithm.cpp    
    DDSensor.cpp    
    DDOptions.cpp
    DDAncillary.cpp
    resam_viirs/RsViirs.cpp
    darktarget/DtLutNetcdf.cpp 
    darktarget/DtMask.cpp    
    darktarget/DtAlgorithm.cpp
    darktarget/DtAlgLand.cpp
    darktarget/DtAlgOcean.cpp    
    deepblue/DbLutNetcdf.cpp 
    deepblue/DbMask.cpp   
    deepblue/DbAlgorithm.cpp   
    deepblue/DbAlgLand.cpp   
    deepblue/DbAlgOcean.cpp   
    deepblue/deepblue_initialize.f90
    deepblue/deepblue_config.f95
    deepblue/calendars.f95 
    deepblue/viirs_ler_luts_nc4.f95
    deepblue/landcover_nc4.f95
    deepblue/get_lut_igbp_land_cover_nc4.f90
    deepblue/get_lut_211sfc_nc4.f90
    deepblue/modis_surface_nc4.f95
    deepblue/viirs_aerosol_luts_nc4.f95
    deepblue/find_v_viirs_nc4.f       
    deepblue/find_v_vegset_nc4.f       
    deepblue/sortrx.f       
    deepblue/pack_412.f       
)

add_executable(generate_dt_luts
    darktarget/generate_dt_luts.cpp
    darktarget/DtLutNetcdf.cpp 
    DDOptions.cpp   
)

add_executable(generate_db_luts
    deepblue/generate_db_luts.cpp
    deepblue/DbLutNetcdf.cpp 
    DDOptions.cpp   
)

add_executable(resam_viirs
    resam_viirs/resam_viirs.cpp
    resam_viirs/RsViirs.cpp 
    DDOptions.cpp   
)

target_link_libraries(dtdb ${LIBS})
target_link_libraries(generate_dt_luts ${LIBS})
target_link_libraries(generate_db_luts ${LIBS})
target_link_libraries(resam_viirs ${LIBS})

install (TARGETS dtdb DESTINATION bin)
install (TARGETS generate_dt_luts DESTINATION bin)
install (TARGETS generate_db_luts DESTINATION bin)
install (TARGETS resam_viirs DESTINATION bin)
