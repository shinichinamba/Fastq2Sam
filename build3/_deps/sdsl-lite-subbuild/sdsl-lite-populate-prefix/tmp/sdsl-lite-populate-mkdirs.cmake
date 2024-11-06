# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/sdsl-lite-src"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/sdsl-lite-build"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/sdsl-lite-subbuild/sdsl-lite-populate-prefix"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/sdsl-lite-subbuild/sdsl-lite-populate-prefix/tmp"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/sdsl-lite-subbuild/sdsl-lite-populate-prefix/src/sdsl-lite-populate-stamp"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/sdsl-lite-subbuild/sdsl-lite-populate-prefix/src"
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/sdsl-lite-subbuild/sdsl-lite-populate-prefix/src/sdsl-lite-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/sdsl-lite-subbuild/sdsl-lite-populate-prefix/src/sdsl-lite-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/sdsl-lite-subbuild/sdsl-lite-populate-prefix/src/sdsl-lite-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
