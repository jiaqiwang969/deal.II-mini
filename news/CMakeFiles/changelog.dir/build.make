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

# Utility rule file for changelog.

# Include the progress variables for this target.
include doc/news/CMakeFiles/changelog.dir/progress.make

doc/news/CMakeFiles/changelog: doc/news/changes.h


doc/news/changes.h: doc/news/changes/header
doc/news/changes.h: doc/news/changes/header_incompatibilities
doc/news/changes.h: doc/news/changes/header_major
doc/news/changes.h: doc/news/changes/header_minor
doc/news/changes.h: doc/news/changes/footer
doc/news/changes.h: doc/news/changes/major
doc/news/changes.h: doc/news/changes/minor
doc/news/changes.h: doc/news/changes/incompatibilities
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating changes.h"
	cd /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/news/changes && /usr/bin/cmake -DOUTPUT_FILE=/home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/news/changes.h -P /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/news/changes/create_changes_h.cmake

changelog: doc/news/CMakeFiles/changelog
changelog: doc/news/changes.h
changelog: doc/news/CMakeFiles/changelog.dir/build.make

.PHONY : changelog

# Rule to build all files generated by this target.
doc/news/CMakeFiles/changelog.dir/build: changelog

.PHONY : doc/news/CMakeFiles/changelog.dir/build

doc/news/CMakeFiles/changelog.dir/clean:
	cd /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/news && $(CMAKE_COMMAND) -P CMakeFiles/changelog.dir/cmake_clean.cmake
.PHONY : doc/news/CMakeFiles/changelog.dir/clean

doc/news/CMakeFiles/changelog.dir/depend:
	cd /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/news /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/news /home/wjq/actions-runner/_work/deal.II-mini/deal.II-mini/dealii/doc/news/CMakeFiles/changelog.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/news/CMakeFiles/changelog.dir/depend
