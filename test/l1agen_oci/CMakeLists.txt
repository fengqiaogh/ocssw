enable_testing()

#include(PrintAllVariables)
#print_all_variables()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l1agen_oci-data"
  VERBOSE
  TEST l1agen_oci-data
  SAVE_DIR dataDir)

# add a test that creates the output directory
add_test(NAME "l1agen_oci-output-rm"
  WORKING_DIRECTORY ${dataDir}
  COMMAND rm -rf output)

add_test(NAME "l1agen_oci-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND mkdir -p output)

################################################################################
# OCI L1A tests
################################################################################
set(l1file "PACE_OCI_2020308009.20201103T145425.L1A.nc")
set(l0file "2020_308_009_IFF.odps")
add_test(NAME "l1agen_oci_${l0file}-test"
    WORKING_DIRECTORY ${dataDir}/output
    COMMAND 
       bash -c "l1agen_oci ../${l0file} 1 --noSPW; \
                cd ${dataDir}; \
                nccmp -m -g -d -f -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${l1file} output/${l1file}")
