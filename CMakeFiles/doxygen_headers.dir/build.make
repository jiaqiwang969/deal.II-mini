# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii

# Include any dependencies generated for this target.
include doc/CMakeFiles/doxygen_headers.dir/depend.make

# Include the progress variables for this target.
include doc/CMakeFiles/doxygen_headers.dir/progress.make

# Include the compile flags for this target's objects.
include doc/CMakeFiles/doxygen_headers.dir/flags.make

doxygen_headers: doc/CMakeFiles/doxygen_headers.dir/build.make

.PHONY : doxygen_headers

# Rule to build all files generated by this target.
doc/CMakeFiles/doxygen_headers.dir/build: doxygen_headers

.PHONY : doc/CMakeFiles/doxygen_headers.dir/build

doc/CMakeFiles/doxygen_headers.dir/requires:

.PHONY : doc/CMakeFiles/doxygen_headers.dir/requires

doc/CMakeFiles/doxygen_headers.dir/clean:
	cd /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc && $(CMAKE_COMMAND) -P CMakeFiles/doxygen_headers.dir/cmake_clean.cmake
.PHONY : doc/CMakeFiles/doxygen_headers.dir/clean

doc/CMakeFiles/doxygen_headers.dir/depend:
	cd /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/CMakeFiles/doxygen_headers.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/CMakeFiles/doxygen_headers.dir/depend
