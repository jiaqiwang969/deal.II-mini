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

# Utility rule file for tutorial_step-77.

# Include the progress variables for this target.
include doc/doxygen/tutorial/CMakeFiles/tutorial_step-77.dir/progress.make

doc/doxygen/tutorial/CMakeFiles/tutorial_step-77: doc/doxygen/tutorial/step-77.h
doc/doxygen/tutorial/CMakeFiles/tutorial_step-77: doc/doxygen/tutorial/step-77.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building doxygen input file for tutorial program <step-77>"

doc/doxygen/tutorial/step-77.h: doc/doxygen/scripts/make_step.pl
doc/doxygen/tutorial/step-77.h: doc/doxygen/scripts/intro2toc
doc/doxygen/tutorial/step-77.h: doc/doxygen/scripts/create_anchors
doc/doxygen/tutorial/step-77.h: doc/doxygen/scripts/program2doxygen
doc/doxygen/tutorial/step-77.h: examples/step-77/step-77.cc
doc/doxygen/tutorial/step-77.h: examples/step-77/doc/intro.dox
doc/doxygen/tutorial/step-77.h: examples/step-77/doc/results.dox
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating step-77.h"
	cd /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/tutorial && /usr/bin/perl /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/scripts/make_step.pl step-77 /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii > /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/tutorial/step-77.h

doc/doxygen/tutorial/step-77.cc: doc/doxygen/scripts/program2plain
doc/doxygen/tutorial/step-77.cc: examples/step-77/step-77.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Generating step-77.cc"
	cd /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/tutorial && /usr/bin/perl /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/scripts/program2plain < /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/examples/step-77/step-77.cc > /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/tutorial/step-77.cc

tutorial_step-77: doc/doxygen/tutorial/CMakeFiles/tutorial_step-77
tutorial_step-77: doc/doxygen/tutorial/step-77.h
tutorial_step-77: doc/doxygen/tutorial/step-77.cc
tutorial_step-77: doc/doxygen/tutorial/CMakeFiles/tutorial_step-77.dir/build.make

.PHONY : tutorial_step-77

# Rule to build all files generated by this target.
doc/doxygen/tutorial/CMakeFiles/tutorial_step-77.dir/build: tutorial_step-77

.PHONY : doc/doxygen/tutorial/CMakeFiles/tutorial_step-77.dir/build

doc/doxygen/tutorial/CMakeFiles/tutorial_step-77.dir/clean:
	cd /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/tutorial && $(CMAKE_COMMAND) -P CMakeFiles/tutorial_step-77.dir/cmake_clean.cmake
.PHONY : doc/doxygen/tutorial/CMakeFiles/tutorial_step-77.dir/clean

doc/doxygen/tutorial/CMakeFiles/tutorial_step-77.dir/depend:
	cd /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/tutorial /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/tutorial /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/doxygen/tutorial/CMakeFiles/tutorial_step-77.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/doxygen/tutorial/CMakeFiles/tutorial_step-77.dir/depend

