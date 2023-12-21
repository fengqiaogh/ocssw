enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l3bin-data"
  VERBOSE
  TEST l3bin-data
  SAVE_DIR dataDir)

# add a test that creates the output directory
add_test(NAME "l3bin-output-rm"
  WORKING_DIRECTORY ${dataDir}
  COMMAND rm -rf output)

add_test(NAME "l3bin-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND mkdir -p output)

# get a list of the netCDF par files in this directory
file(GLOB files "${dataDir}/*.nc.par")

# loop through all netCDF par files found and make tests
foreach(filename ${files})

  GET_FILENAME_COMPONENT(parfile ${filename} NAME)
  STRING(REGEX REPLACE ".par" "" outfile ${parfile})

  add_test(NAME "l3bin-${outfile}-test"
    WORKING_DIRECTORY ${dataDir}
    COMMAND 
      bash -c "l3bin par=${parfile}; \
               nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} ${outfile} output/${outfile}; \
               compare_l3b ${outfile} output/${outfile}")

endforeach(filename)
