# CMake generated Testfile for 
# Source directory: /mnt/c/Users/sinnh/Dropbox/Fastq2Sam
# Build directory: /mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_run "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build/fastq2sam" "-1" "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/test/1.fq" "-2" "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/test/2.fq" "-o" "test.sam" "-b" "2" "-n" "sample1" "--no-pg")
set_tests_properties(test_run PROPERTIES  WORKING_DIRECTORY "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build" _BACKTRACE_TRIPLES "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/CMakeLists.txt;35;add_test;/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/CMakeLists.txt;0;")
add_test(diff_after_test_run "diff" "-q" "test.sam" "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/test/out.noPG.sam")
set_tests_properties(diff_after_test_run PROPERTIES  WORKING_DIRECTORY "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build" _BACKTRACE_TRIPLES "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/CMakeLists.txt;40;add_test;/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/CMakeLists.txt;0;")
add_test(clean "rm" "-f" "test.sam")
set_tests_properties(clean PROPERTIES  WORKING_DIRECTORY "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build" _BACKTRACE_TRIPLES "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/CMakeLists.txt;45;add_test;/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/CMakeLists.txt;0;")
subdirs("_deps/yaml-cpp_fetch_content-build")
