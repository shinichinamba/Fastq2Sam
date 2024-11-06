# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if(EXISTS "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/cereal-populate-gitinfo.txt" AND EXISTS "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps" AND
  "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/cereal-populate-gitinfo.txt" IS_NEWER_THAN "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps")
  message(STATUS
    "Avoiding repeated git clone, stamp file is up to date: "
    "'/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/cereal-populate-gitinfo.txt'"
  )
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-src"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"
            clone --no-checkout --origin "v1.3.2" -c http.sslVerify=true "https://github.com/USCiLab/cereal.git" "advice.detachedHead=false"
    WORKING_DIRECTORY "cereal-src"
    RESULT_VARIABLE error_code
  )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once: ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/USCiLab/cereal.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"
          checkout "YES" --
  WORKING_DIRECTORY "cereal-src/advice.detachedHead=false"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'YES'")
endif()

set(init_submodules origin)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git" -c;http.sslVerify=true
            submodule update TRUE --init --recursive
    WORKING_DIRECTORY "cereal-src/advice.detachedHead=false"
    RESULT_VARIABLE error_code
  )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: 'cereal-src/advice.detachedHead=false'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps" "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/cereal-populate-gitinfo.txt"
  RESULT_VARIABLE error_code
)
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/cereal-subbuild/cereal-populate-prefix/src/cereal-populate-stamp/cereal-populate-gitinfo.txt'")
endif()
