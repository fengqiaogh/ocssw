enable_testing()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/geogen_modis-data"
  VERBOSE
  TEST geogen_modis-data
  SAVE_DIR dataDir)

# add a test that creates the output directory
add_test(NAME "geogen_modis-output-rm"
  WORKING_DIRECTORY ${dataDir}
  COMMAND rm -rf output)

add_test(NAME "geogen_modis-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND mkdir -p output)

################################################################################
# loop through L1A files making GEO files
################################################################################
file(GLOB files "${dataDir}/*L1A_LAC*")
foreach(filename ${files})

  GET_FILENAME_COMPONENT(l1afile ${filename} NAME)
  STRING(REGEX REPLACE "L1A_LAC" "GEO" geofile ${l1afile})
  STRING(CONCAT pcf ${dataDir}/${geofile} ".pcf")

  add_test(NAME "geogen_modis-${geofile}-test"
    WORKING_DIRECTORY ${dataDir}
    COMMAND 
      bash -c "cp ${l1afile} output; \
               cd ${dataDir}/output; \
               python ../geogen_modis.py ${pcf}; \
               cd ${dataDir}; \
               compare_hdf4 -G \"CoreMetadata.0,utcpole File Header\" ${geofile} output/${geofile}")

endforeach(filename)

################################################################################
# loop through L1A files making L1B files
################################################################################
foreach(filename ${files})

  GET_FILENAME_COMPONENT(l1afile ${filename} NAME)
  STRING(REGEX REPLACE "L1A_LAC" "L1B_LAC" okmfile ${l1afile})
  STRING(REGEX REPLACE "L1A_LAC" "L1B_HKM" hkmfile ${l1afile})
  STRING(REGEX REPLACE "L1A_LAC" "L1B_QKM" qkmfile ${l1afile})
  STRING(REGEX REPLACE "L1A_LAC" "L1B_OBC" obcfile ${l1afile})
  STRING(REGEX REPLACE "L1A_LAC" "GEO" geofile ${l1afile})
  STRING(REGEX REPLACE "L1A_LAC." "" testfile ${l1afile})
  if(l1afile MATCHES "^A")
    set(lutdir "${dataDir}/luts/aqua")
  else()
    set(lutdir "${dataDir}/luts/terra")
  endif()

  add_test(NAME "l1bgen_modis-${testfile}-test"
    WORKING_DIRECTORY ${dataDir}/output
    COMMAND 
      bash -c "modis_L1B --keep-obc --lutdir=${lutdir} --okm ${okmfile} --hkm ${hkmfile} --qkm ${qkmfile} --obc ${obcfile} ${l1afile} ../${geofile}; \
               cd ${dataDir}; \
               compare_hdf4 -G \"CoreMetadata.0,ArchiveMetadata.0\" ${okmfile} output/${okmfile};\
               compare_hdf4 -G \"CoreMetadata.0,ArchiveMetadata.0\" ${hkmfile} output/${hkmfile};\
               compare_hdf4 -G \"CoreMetadata.0,ArchiveMetadata.0\" ${qkmfile} output/${qkmfile};\
               compare_hdf4 -G \"CoreMetadata.0,ArchiveMetadata.0\" ${obcfile} output/${obcfile}")

endforeach(filename)
