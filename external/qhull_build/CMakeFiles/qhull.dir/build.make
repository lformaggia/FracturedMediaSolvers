# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/cmake-gui

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/forma/workspace2/eni/eniReservoir/external/qhull-2012.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/forma/workspace2/eni/eniReservoir/external/qhull_build

# Include any dependencies generated for this target.
include CMakeFiles/qhull.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/qhull.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/qhull.dir/flags.make

CMakeFiles/qhull.dir/src/qhull/unix.c.o: CMakeFiles/qhull.dir/flags.make
CMakeFiles/qhull.dir/src/qhull/unix.c.o: /home/forma/workspace2/eni/eniReservoir/external/qhull-2012.1/src/qhull/unix.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/forma/workspace2/eni/eniReservoir/external/qhull_build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/qhull.dir/src/qhull/unix.c.o"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/qhull.dir/src/qhull/unix.c.o   -c /home/forma/workspace2/eni/eniReservoir/external/qhull-2012.1/src/qhull/unix.c

CMakeFiles/qhull.dir/src/qhull/unix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/qhull.dir/src/qhull/unix.c.i"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/forma/workspace2/eni/eniReservoir/external/qhull-2012.1/src/qhull/unix.c > CMakeFiles/qhull.dir/src/qhull/unix.c.i

CMakeFiles/qhull.dir/src/qhull/unix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/qhull.dir/src/qhull/unix.c.s"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/forma/workspace2/eni/eniReservoir/external/qhull-2012.1/src/qhull/unix.c -o CMakeFiles/qhull.dir/src/qhull/unix.c.s

CMakeFiles/qhull.dir/src/qhull/unix.c.o.requires:
.PHONY : CMakeFiles/qhull.dir/src/qhull/unix.c.o.requires

CMakeFiles/qhull.dir/src/qhull/unix.c.o.provides: CMakeFiles/qhull.dir/src/qhull/unix.c.o.requires
	$(MAKE) -f CMakeFiles/qhull.dir/build.make CMakeFiles/qhull.dir/src/qhull/unix.c.o.provides.build
.PHONY : CMakeFiles/qhull.dir/src/qhull/unix.c.o.provides

CMakeFiles/qhull.dir/src/qhull/unix.c.o.provides.build: CMakeFiles/qhull.dir/src/qhull/unix.c.o

# Object files for target qhull
qhull_OBJECTS = \
"CMakeFiles/qhull.dir/src/qhull/unix.c.o"

# External object files for target qhull
qhull_EXTERNAL_OBJECTS =

qhull: CMakeFiles/qhull.dir/src/qhull/unix.c.o
qhull: libqhullstatic.a
qhull: CMakeFiles/qhull.dir/build.make
qhull: CMakeFiles/qhull.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable qhull"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/qhull.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/qhull.dir/build: qhull
.PHONY : CMakeFiles/qhull.dir/build

CMakeFiles/qhull.dir/requires: CMakeFiles/qhull.dir/src/qhull/unix.c.o.requires
.PHONY : CMakeFiles/qhull.dir/requires

CMakeFiles/qhull.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/qhull.dir/cmake_clean.cmake
.PHONY : CMakeFiles/qhull.dir/clean

CMakeFiles/qhull.dir/depend:
	cd /home/forma/workspace2/eni/eniReservoir/external/qhull_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/forma/workspace2/eni/eniReservoir/external/qhull-2012.1 /home/forma/workspace2/eni/eniReservoir/external/qhull-2012.1 /home/forma/workspace2/eni/eniReservoir/external/qhull_build /home/forma/workspace2/eni/eniReservoir/external/qhull_build /home/forma/workspace2/eni/eniReservoir/external/qhull_build/CMakeFiles/qhull.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/qhull.dir/depend

