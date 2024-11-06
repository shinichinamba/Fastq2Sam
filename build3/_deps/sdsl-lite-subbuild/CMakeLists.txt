# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.27.0-rc3)

# Reject any attempt to use a toolchain file. We must not use one because
# we could be downloading it here. If the CMAKE_TOOLCHAIN_FILE environment
# variable is set, the cache variable will have been initialized from it.
unset(CMAKE_TOOLCHAIN_FILE CACHE)
unset(ENV{CMAKE_TOOLCHAIN_FILE})

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(sdsl-lite-populate NONE)



include(ExternalProject)
ExternalProject_Add(sdsl-lite-populate
                     "UPDATE_DISCONNECTED" "False" "DOWNLOAD_EXTRACT_TIMESTAMP" "FALSE" "YES" "URL" "https://github.com/xxsds/sdsl-lite/archive/14cd017027ea742353fc5b500d1cb1d95896b77e.tar.gz"
                    SOURCE_DIR          "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/sdsl-lite-src"
                    BINARY_DIR          "/mnt/c/Users/sinnh/Dropbox/Fastq2Sam/build3/_deps/sdsl-lite-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
                    USES_TERMINAL_PATCH     YES
)


