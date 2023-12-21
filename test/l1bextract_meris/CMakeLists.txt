enable_testing()

cdata(READ data_files "${CMAKE_CURRENT_LIST_DIR}/CDataList.txt"
  GENERATE "${CMAKE_CURRENT_BINARY_DIR}/l1bextract_meris-data"
  VERBOSE
  TEST l1bextract_meris-data
  SAVE_DIR dataDir)

# add a test that creates the output directory
add_test(NAME "l1bextract_meris-output-rm"
  WORKING_DIRECTORY ${dataDir}
  COMMAND rm -rf output)

add_test(NAME "l1bextract_meris-output"
  WORKING_DIRECTORY ${dataDir}
  COMMAND mkdir -p output)

add_test(NAME "l1bextract_meris-test"
  WORKING_DIRECTORY ${dataDir}
  COMMAND 
    bash -c "l1bextract_meris MER_RR__1PNACR20060921_150933_000026302051_00240_23841_0000.N1 400 900 5377 5900 output/extract.N1; \
             l2gen par=extract.par; \
             nccmp -m -g -d -f -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} extract.L2 output/extract.L2")

# test the SAFE extractor
add_test(NAME "l1bextract_meris_safe-test"
  WORKING_DIRECTORY ${dataDir}
  COMMAND 
    bash -c "l1bextract_safe_nc --spixl 300 --epixl 900 --sline 6000 --eline 6200 ENV_ME_1_RRG____20061115T162732_20061115T171127_________________2636_053_026______DSI_R_NT____.SEN3 output/safe.extract; \
             cd ${dataDir}/output/safe.extract; \
             l2gen atmocor=off l2prod=Lt_490 ifile=xfdumanifest.xml ofile=meris.safe.extract.L2.nc; \

             cd ${dataDir}; \
             nccmp -m -g -d -f -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} meris.safe.extract.L2.nc output/safe.extract/meris.safe.extract.L2.nc")
