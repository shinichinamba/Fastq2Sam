# CMake generated Testfile for 
# Source directory: /Users/snamba/Dropbox/Fastq2Sam
# Build directory: /Users/snamba/Dropbox/Fastq2Sam/build_mac
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_run "/Users/snamba/Dropbox/Fastq2Sam/build_mac/fastq2sam" "-1" "/Users/snamba/Dropbox/Fastq2Sam/test/1.fq" "-2" "/Users/snamba/Dropbox/Fastq2Sam/test/2.fq" "-o" "test.sam" "-b" "2" "-n" "sample1" "--no-pg")
set_tests_properties(test_run PROPERTIES  WORKING_DIRECTORY "/Users/snamba/Dropbox/Fastq2Sam/build_mac" _BACKTRACE_TRIPLES "/Users/snamba/Dropbox/Fastq2Sam/CMakeLists.txt;29;add_test;/Users/snamba/Dropbox/Fastq2Sam/CMakeLists.txt;0;")
add_test(diff_after_test_run "diff" "-q" "test.sam" "/Users/snamba/Dropbox/Fastq2Sam/test/out.noPG.sam")
set_tests_properties(diff_after_test_run PROPERTIES  WORKING_DIRECTORY "/Users/snamba/Dropbox/Fastq2Sam/build_mac" _BACKTRACE_TRIPLES "/Users/snamba/Dropbox/Fastq2Sam/CMakeLists.txt;34;add_test;/Users/snamba/Dropbox/Fastq2Sam/CMakeLists.txt;0;")
add_test(clean "rm" "-f" "test.sam")
set_tests_properties(clean PROPERTIES  WORKING_DIRECTORY "/Users/snamba/Dropbox/Fastq2Sam/build_mac" _BACKTRACE_TRIPLES "/Users/snamba/Dropbox/Fastq2Sam/CMakeLists.txt;39;add_test;/Users/snamba/Dropbox/Fastq2Sam/CMakeLists.txt;0;")
subdirs("_deps/yaml-cpp_fetch_content-build")
