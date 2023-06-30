cmake_minimum_required (VERSION 3.4)
project (fastq2sam CXX)
set (CMAKE_BUILD_TYPE Release)

# set (BUILD_SHARED_LIBS false)
# set (CMAKE_EXE_LINKER_FLAGS -static)
# set (ZLIB_USE_STATIC_LIBS ON)
# set (BZIP2_STATIC_LIBRARIES ON)

# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/seqan3/build_system")
# add the Sharg Parser to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/sharg-parser/build_system")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)
find_package (sharg 1.0 REQUIRED)

# build app with seqan3
add_executable (fastq2sam src/fastq2sam.cpp) # src/fastq2sam.seqan3.2.cpp
target_link_libraries (fastq2sam seqan3::seqan3 sharg::sharg) #  -static
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    target_compile_options(fastq2sam PUBLIC -O3 -Wall)
endif()
target_compile_features(fastq2sam PUBLIC cxx_std_20)

# ctest
enable_testing()
add_test(
    NAME test_run
    COMMAND $<TARGET_FILE:fastq2sam> -1 "${CMAKE_CURRENT_SOURCE_DIR}/test/1.fq" -2 "${CMAKE_CURRENT_SOURCE_DIR}/test/2.fq" -o test.sam -b 2 -n sample1 --no-pg
    WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
)
add_test(
    NAME diff_after_test_run
    COMMAND diff -q test.sam "${CMAKE_CURRENT_SOURCE_DIR}/test/out.noPG.sam"
    WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
)
add_test(
    NAME clean
    COMMAND rm -f test.sam
    WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
)
