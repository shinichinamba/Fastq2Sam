# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.17.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.17.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/snamba/Dropbox/Fastq2Sam

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/snamba/Dropbox/Fastq2Sam/build_mac

# Include any dependencies generated for this target.
include CMakeFiles/fastq2sam.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/fastq2sam.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fastq2sam.dir/flags.make

CMakeFiles/fastq2sam.dir/src/bamhash.cpp.o: CMakeFiles/fastq2sam.dir/flags.make
CMakeFiles/fastq2sam.dir/src/bamhash.cpp.o: ../src/bamhash.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/snamba/Dropbox/Fastq2Sam/build_mac/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/fastq2sam.dir/src/bamhash.cpp.o"
	/usr/local/bin/g++-13  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fastq2sam.dir/src/bamhash.cpp.o -c /Users/snamba/Dropbox/Fastq2Sam/src/bamhash.cpp

CMakeFiles/fastq2sam.dir/src/bamhash.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fastq2sam.dir/src/bamhash.cpp.i"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/snamba/Dropbox/Fastq2Sam/src/bamhash.cpp > CMakeFiles/fastq2sam.dir/src/bamhash.cpp.i

CMakeFiles/fastq2sam.dir/src/bamhash.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fastq2sam.dir/src/bamhash.cpp.s"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/snamba/Dropbox/Fastq2Sam/src/bamhash.cpp -o CMakeFiles/fastq2sam.dir/src/bamhash.cpp.s

CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.o: CMakeFiles/fastq2sam.dir/flags.make
CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.o: ../src/convert_phred_score.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/snamba/Dropbox/Fastq2Sam/build_mac/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.o"
	/usr/local/bin/g++-13  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.o -c /Users/snamba/Dropbox/Fastq2Sam/src/convert_phred_score.cpp

CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.i"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/snamba/Dropbox/Fastq2Sam/src/convert_phred_score.cpp > CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.i

CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.s"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/snamba/Dropbox/Fastq2Sam/src/convert_phred_score.cpp -o CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.s

CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.o: CMakeFiles/fastq2sam.dir/flags.make
CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.o: ../src/fastq_metadata.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/snamba/Dropbox/Fastq2Sam/build_mac/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.o"
	/usr/local/bin/g++-13  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.o -c /Users/snamba/Dropbox/Fastq2Sam/src/fastq_metadata.cpp

CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.i"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/snamba/Dropbox/Fastq2Sam/src/fastq_metadata.cpp > CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.i

CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.s"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/snamba/Dropbox/Fastq2Sam/src/fastq_metadata.cpp -o CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.s

CMakeFiles/fastq2sam.dir/src/id_parser.cpp.o: CMakeFiles/fastq2sam.dir/flags.make
CMakeFiles/fastq2sam.dir/src/id_parser.cpp.o: ../src/id_parser.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/snamba/Dropbox/Fastq2Sam/build_mac/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/fastq2sam.dir/src/id_parser.cpp.o"
	/usr/local/bin/g++-13  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fastq2sam.dir/src/id_parser.cpp.o -c /Users/snamba/Dropbox/Fastq2Sam/src/id_parser.cpp

CMakeFiles/fastq2sam.dir/src/id_parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fastq2sam.dir/src/id_parser.cpp.i"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/snamba/Dropbox/Fastq2Sam/src/id_parser.cpp > CMakeFiles/fastq2sam.dir/src/id_parser.cpp.i

CMakeFiles/fastq2sam.dir/src/id_parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fastq2sam.dir/src/id_parser.cpp.s"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/snamba/Dropbox/Fastq2Sam/src/id_parser.cpp -o CMakeFiles/fastq2sam.dir/src/id_parser.cpp.s

CMakeFiles/fastq2sam.dir/src/phred.cpp.o: CMakeFiles/fastq2sam.dir/flags.make
CMakeFiles/fastq2sam.dir/src/phred.cpp.o: ../src/phred.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/snamba/Dropbox/Fastq2Sam/build_mac/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/fastq2sam.dir/src/phred.cpp.o"
	/usr/local/bin/g++-13  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fastq2sam.dir/src/phred.cpp.o -c /Users/snamba/Dropbox/Fastq2Sam/src/phred.cpp

CMakeFiles/fastq2sam.dir/src/phred.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fastq2sam.dir/src/phred.cpp.i"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/snamba/Dropbox/Fastq2Sam/src/phred.cpp > CMakeFiles/fastq2sam.dir/src/phred.cpp.i

CMakeFiles/fastq2sam.dir/src/phred.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fastq2sam.dir/src/phred.cpp.s"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/snamba/Dropbox/Fastq2Sam/src/phred.cpp -o CMakeFiles/fastq2sam.dir/src/phred.cpp.s

CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.o: CMakeFiles/fastq2sam.dir/flags.make
CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.o: ../src/scan_fastq.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/snamba/Dropbox/Fastq2Sam/build_mac/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.o"
	/usr/local/bin/g++-13  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.o -c /Users/snamba/Dropbox/Fastq2Sam/src/scan_fastq.cpp

CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.i"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/snamba/Dropbox/Fastq2Sam/src/scan_fastq.cpp > CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.i

CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.s"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/snamba/Dropbox/Fastq2Sam/src/scan_fastq.cpp -o CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.s

CMakeFiles/fastq2sam.dir/src/split.cpp.o: CMakeFiles/fastq2sam.dir/flags.make
CMakeFiles/fastq2sam.dir/src/split.cpp.o: ../src/split.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/snamba/Dropbox/Fastq2Sam/build_mac/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/fastq2sam.dir/src/split.cpp.o"
	/usr/local/bin/g++-13  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fastq2sam.dir/src/split.cpp.o -c /Users/snamba/Dropbox/Fastq2Sam/src/split.cpp

CMakeFiles/fastq2sam.dir/src/split.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fastq2sam.dir/src/split.cpp.i"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/snamba/Dropbox/Fastq2Sam/src/split.cpp > CMakeFiles/fastq2sam.dir/src/split.cpp.i

CMakeFiles/fastq2sam.dir/src/split.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fastq2sam.dir/src/split.cpp.s"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/snamba/Dropbox/Fastq2Sam/src/split.cpp -o CMakeFiles/fastq2sam.dir/src/split.cpp.s

CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.o: CMakeFiles/fastq2sam.dir/flags.make
CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.o: ../src/fastq2sam.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/snamba/Dropbox/Fastq2Sam/build_mac/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.o"
	/usr/local/bin/g++-13  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.o -c /Users/snamba/Dropbox/Fastq2Sam/src/fastq2sam.cpp

CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.i"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/snamba/Dropbox/Fastq2Sam/src/fastq2sam.cpp > CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.i

CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.s"
	/usr/local/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/snamba/Dropbox/Fastq2Sam/src/fastq2sam.cpp -o CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.s

# Object files for target fastq2sam
fastq2sam_OBJECTS = \
"CMakeFiles/fastq2sam.dir/src/bamhash.cpp.o" \
"CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.o" \
"CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.o" \
"CMakeFiles/fastq2sam.dir/src/id_parser.cpp.o" \
"CMakeFiles/fastq2sam.dir/src/phred.cpp.o" \
"CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.o" \
"CMakeFiles/fastq2sam.dir/src/split.cpp.o" \
"CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.o"

# External object files for target fastq2sam
fastq2sam_EXTERNAL_OBJECTS =

fastq2sam: CMakeFiles/fastq2sam.dir/src/bamhash.cpp.o
fastq2sam: CMakeFiles/fastq2sam.dir/src/convert_phred_score.cpp.o
fastq2sam: CMakeFiles/fastq2sam.dir/src/fastq_metadata.cpp.o
fastq2sam: CMakeFiles/fastq2sam.dir/src/id_parser.cpp.o
fastq2sam: CMakeFiles/fastq2sam.dir/src/phred.cpp.o
fastq2sam: CMakeFiles/fastq2sam.dir/src/scan_fastq.cpp.o
fastq2sam: CMakeFiles/fastq2sam.dir/src/split.cpp.o
fastq2sam: CMakeFiles/fastq2sam.dir/src/fastq2sam.cpp.o
fastq2sam: CMakeFiles/fastq2sam.dir/build.make
fastq2sam: /usr/local/Cellar/openssl@3/3.1.1_1/lib/libcrypto.dylib
fastq2sam: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.0.sdk/usr/lib/libz.tbd
fastq2sam: _deps/yaml-cpp_fetch_content-build/libyaml-cppd.a
fastq2sam: CMakeFiles/fastq2sam.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/snamba/Dropbox/Fastq2Sam/build_mac/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX executable fastq2sam"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fastq2sam.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fastq2sam.dir/build: fastq2sam

.PHONY : CMakeFiles/fastq2sam.dir/build

CMakeFiles/fastq2sam.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fastq2sam.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fastq2sam.dir/clean

CMakeFiles/fastq2sam.dir/depend:
	cd /Users/snamba/Dropbox/Fastq2Sam/build_mac && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/snamba/Dropbox/Fastq2Sam /Users/snamba/Dropbox/Fastq2Sam /Users/snamba/Dropbox/Fastq2Sam/build_mac /Users/snamba/Dropbox/Fastq2Sam/build_mac /Users/snamba/Dropbox/Fastq2Sam/build_mac/CMakeFiles/fastq2sam.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fastq2sam.dir/depend

