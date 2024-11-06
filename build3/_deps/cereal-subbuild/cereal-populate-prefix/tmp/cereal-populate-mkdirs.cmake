# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-src"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-build"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix/tmp"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix/src"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
