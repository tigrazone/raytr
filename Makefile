# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_COMMAND = C:/CMake/bin/cmake.exe

# The command to remove a file.
RM = C:/CMake/bin/cmake.exe -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:/msys64/home/thr/raytr-master

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:/msys64/home/thr/raytr-master

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	C:/CMake/bin/cmake-gui.exe -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	C:/CMake/bin/cmake.exe -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start C:/msys64/home/thr/raytr-master/CMakeFiles C:/msys64/home/thr/raytr-master/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start C:/msys64/home/thr/raytr-master/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named raytr

# Build rule for target.
raytr: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 raytr
.PHONY : raytr

# fast build rule for target.
raytr/fast:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/build
.PHONY : raytr/fast

bbox_tree.obj: bbox_tree.cpp.obj

.PHONY : bbox_tree.obj

# target to build an object file
bbox_tree.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/bbox_tree.cpp.obj
.PHONY : bbox_tree.cpp.obj

bbox_tree.i: bbox_tree.cpp.i

.PHONY : bbox_tree.i

# target to preprocess a source file
bbox_tree.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/bbox_tree.cpp.i
.PHONY : bbox_tree.cpp.i

bbox_tree.s: bbox_tree.cpp.s

.PHONY : bbox_tree.s

# target to generate assembly for a file
bbox_tree.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/bbox_tree.cpp.s
.PHONY : bbox_tree.cpp.s

ext/SimpleObject.obj: ext/SimpleObject.cpp.obj

.PHONY : ext/SimpleObject.obj

# target to build an object file
ext/SimpleObject.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/ext/SimpleObject.cpp.obj
.PHONY : ext/SimpleObject.cpp.obj

ext/SimpleObject.i: ext/SimpleObject.cpp.i

.PHONY : ext/SimpleObject.i

# target to preprocess a source file
ext/SimpleObject.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/ext/SimpleObject.cpp.i
.PHONY : ext/SimpleObject.cpp.i

ext/SimpleObject.s: ext/SimpleObject.cpp.s

.PHONY : ext/SimpleObject.s

# target to generate assembly for a file
ext/SimpleObject.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/ext/SimpleObject.cpp.s
.PHONY : ext/SimpleObject.cpp.s

ext/dSFMT/dSFMT.obj: ext/dSFMT/dSFMT.cpp.obj

.PHONY : ext/dSFMT/dSFMT.obj

# target to build an object file
ext/dSFMT/dSFMT.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/ext/dSFMT/dSFMT.cpp.obj
.PHONY : ext/dSFMT/dSFMT.cpp.obj

ext/dSFMT/dSFMT.i: ext/dSFMT/dSFMT.cpp.i

.PHONY : ext/dSFMT/dSFMT.i

# target to preprocess a source file
ext/dSFMT/dSFMT.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/ext/dSFMT/dSFMT.cpp.i
.PHONY : ext/dSFMT/dSFMT.cpp.i

ext/dSFMT/dSFMT.s: ext/dSFMT/dSFMT.cpp.s

.PHONY : ext/dSFMT/dSFMT.s

# target to generate assembly for a file
ext/dSFMT/dSFMT.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/ext/dSFMT/dSFMT.cpp.s
.PHONY : ext/dSFMT/dSFMT.cpp.s

ext/tiny_obj_loader.obj: ext/tiny_obj_loader.cc.obj

.PHONY : ext/tiny_obj_loader.obj

# target to build an object file
ext/tiny_obj_loader.cc.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/ext/tiny_obj_loader.cc.obj
.PHONY : ext/tiny_obj_loader.cc.obj

ext/tiny_obj_loader.i: ext/tiny_obj_loader.cc.i

.PHONY : ext/tiny_obj_loader.i

# target to preprocess a source file
ext/tiny_obj_loader.cc.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/ext/tiny_obj_loader.cc.i
.PHONY : ext/tiny_obj_loader.cc.i

ext/tiny_obj_loader.s: ext/tiny_obj_loader.cc.s

.PHONY : ext/tiny_obj_loader.s

# target to generate assembly for a file
ext/tiny_obj_loader.cc.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/ext/tiny_obj_loader.cc.s
.PHONY : ext/tiny_obj_loader.cc.s

image.obj: image.cpp.obj

.PHONY : image.obj

# target to build an object file
image.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/image.cpp.obj
.PHONY : image.cpp.obj

image.i: image.cpp.i

.PHONY : image.i

# target to preprocess a source file
image.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/image.cpp.i
.PHONY : image.cpp.i

image.s: image.cpp.s

.PHONY : image.s

# target to generate assembly for a file
image.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/image.cpp.s
.PHONY : image.cpp.s

main.obj: main.cpp.obj

.PHONY : main.obj

# target to build an object file
main.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/main.cpp.obj
.PHONY : main.cpp.obj

main.i: main.cpp.i

.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s

.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/main.cpp.s
.PHONY : main.cpp.s

meshsimp1.obj: meshsimp1.cpp.obj

.PHONY : meshsimp1.obj

# target to build an object file
meshsimp1.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/meshsimp1.cpp.obj
.PHONY : meshsimp1.cpp.obj

meshsimp1.i: meshsimp1.cpp.i

.PHONY : meshsimp1.i

# target to preprocess a source file
meshsimp1.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/meshsimp1.cpp.i
.PHONY : meshsimp1.cpp.i

meshsimp1.s: meshsimp1.cpp.s

.PHONY : meshsimp1.s

# target to generate assembly for a file
meshsimp1.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/meshsimp1.cpp.s
.PHONY : meshsimp1.cpp.s

objekt.obj: objekt.cpp.obj

.PHONY : objekt.obj

# target to build an object file
objekt.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/objekt.cpp.obj
.PHONY : objekt.cpp.obj

objekt.i: objekt.cpp.i

.PHONY : objekt.i

# target to preprocess a source file
objekt.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/objekt.cpp.i
.PHONY : objekt.cpp.i

objekt.s: objekt.cpp.s

.PHONY : objekt.s

# target to generate assembly for a file
objekt.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/objekt.cpp.s
.PHONY : objekt.cpp.s

objloader.obj: objloader.cpp.obj

.PHONY : objloader.obj

# target to build an object file
objloader.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/objloader.cpp.obj
.PHONY : objloader.cpp.obj

objloader.i: objloader.cpp.i

.PHONY : objloader.i

# target to preprocess a source file
objloader.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/objloader.cpp.i
.PHONY : objloader.cpp.i

objloader.s: objloader.cpp.s

.PHONY : objloader.s

# target to generate assembly for a file
objloader.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/objloader.cpp.s
.PHONY : objloader.cpp.s

pathtracer.obj: pathtracer.cpp.obj

.PHONY : pathtracer.obj

# target to build an object file
pathtracer.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/pathtracer.cpp.obj
.PHONY : pathtracer.cpp.obj

pathtracer.i: pathtracer.cpp.i

.PHONY : pathtracer.i

# target to preprocess a source file
pathtracer.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/pathtracer.cpp.i
.PHONY : pathtracer.cpp.i

pathtracer.s: pathtracer.cpp.s

.PHONY : pathtracer.s

# target to generate assembly for a file
pathtracer.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/pathtracer.cpp.s
.PHONY : pathtracer.cpp.s

photon_tree.obj: photon_tree.cpp.obj

.PHONY : photon_tree.obj

# target to build an object file
photon_tree.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/photon_tree.cpp.obj
.PHONY : photon_tree.cpp.obj

photon_tree.i: photon_tree.cpp.i

.PHONY : photon_tree.i

# target to preprocess a source file
photon_tree.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/photon_tree.cpp.i
.PHONY : photon_tree.cpp.i

photon_tree.s: photon_tree.cpp.s

.PHONY : photon_tree.s

# target to generate assembly for a file
photon_tree.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/photon_tree.cpp.s
.PHONY : photon_tree.cpp.s

pm.obj: pm.cpp.obj

.PHONY : pm.obj

# target to build an object file
pm.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/pm.cpp.obj
.PHONY : pm.cpp.obj

pm.i: pm.cpp.i

.PHONY : pm.i

# target to preprocess a source file
pm.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/pm.cpp.i
.PHONY : pm.cpp.i

pm.s: pm.cpp.s

.PHONY : pm.s

# target to generate assembly for a file
pm.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/pm.cpp.s
.PHONY : pm.cpp.s

primitives.obj: primitives.cpp.obj

.PHONY : primitives.obj

# target to build an object file
primitives.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/primitives.cpp.obj
.PHONY : primitives.cpp.obj

primitives.i: primitives.cpp.i

.PHONY : primitives.i

# target to preprocess a source file
primitives.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/primitives.cpp.i
.PHONY : primitives.cpp.i

primitives.s: primitives.cpp.s

.PHONY : primitives.s

# target to generate assembly for a file
primitives.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/primitives.cpp.s
.PHONY : primitives.cpp.s

random.obj: random.cpp.obj

.PHONY : random.obj

# target to build an object file
random.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/random.cpp.obj
.PHONY : random.cpp.obj

random.i: random.cpp.i

.PHONY : random.i

# target to preprocess a source file
random.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/random.cpp.i
.PHONY : random.cpp.i

random.s: random.cpp.s

.PHONY : random.s

# target to generate assembly for a file
random.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/random.cpp.s
.PHONY : random.cpp.s

scene.obj: scene.cpp.obj

.PHONY : scene.obj

# target to build an object file
scene.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/scene.cpp.obj
.PHONY : scene.cpp.obj

scene.i: scene.cpp.i

.PHONY : scene.i

# target to preprocess a source file
scene.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/scene.cpp.i
.PHONY : scene.cpp.i

scene.s: scene.cpp.s

.PHONY : scene.s

# target to generate assembly for a file
scene.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/scene.cpp.s
.PHONY : scene.cpp.s

shade.obj: shade.cpp.obj

.PHONY : shade.obj

# target to build an object file
shade.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/shade.cpp.obj
.PHONY : shade.cpp.obj

shade.i: shade.cpp.i

.PHONY : shade.i

# target to preprocess a source file
shade.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/shade.cpp.i
.PHONY : shade.cpp.i

shade.s: shade.cpp.s

.PHONY : shade.s

# target to generate assembly for a file
shade.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/shade.cpp.s
.PHONY : shade.cpp.s

spectrum.obj: spectrum.cpp.obj

.PHONY : spectrum.obj

# target to build an object file
spectrum.cpp.obj:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/spectrum.cpp.obj
.PHONY : spectrum.cpp.obj

spectrum.i: spectrum.cpp.i

.PHONY : spectrum.i

# target to preprocess a source file
spectrum.cpp.i:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/spectrum.cpp.i
.PHONY : spectrum.cpp.i

spectrum.s: spectrum.cpp.s

.PHONY : spectrum.s

# target to generate assembly for a file
spectrum.cpp.s:
	$(MAKE) -f CMakeFiles/raytr.dir/build.make CMakeFiles/raytr.dir/spectrum.cpp.s
.PHONY : spectrum.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... raytr"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... bbox_tree.obj"
	@echo "... bbox_tree.i"
	@echo "... bbox_tree.s"
	@echo "... ext/SimpleObject.obj"
	@echo "... ext/SimpleObject.i"
	@echo "... ext/SimpleObject.s"
	@echo "... ext/dSFMT/dSFMT.obj"
	@echo "... ext/dSFMT/dSFMT.i"
	@echo "... ext/dSFMT/dSFMT.s"
	@echo "... ext/tiny_obj_loader.obj"
	@echo "... ext/tiny_obj_loader.i"
	@echo "... ext/tiny_obj_loader.s"
	@echo "... image.obj"
	@echo "... image.i"
	@echo "... image.s"
	@echo "... main.obj"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... meshsimp1.obj"
	@echo "... meshsimp1.i"
	@echo "... meshsimp1.s"
	@echo "... objekt.obj"
	@echo "... objekt.i"
	@echo "... objekt.s"
	@echo "... objloader.obj"
	@echo "... objloader.i"
	@echo "... objloader.s"
	@echo "... pathtracer.obj"
	@echo "... pathtracer.i"
	@echo "... pathtracer.s"
	@echo "... photon_tree.obj"
	@echo "... photon_tree.i"
	@echo "... photon_tree.s"
	@echo "... pm.obj"
	@echo "... pm.i"
	@echo "... pm.s"
	@echo "... primitives.obj"
	@echo "... primitives.i"
	@echo "... primitives.s"
	@echo "... random.obj"
	@echo "... random.i"
	@echo "... random.s"
	@echo "... scene.obj"
	@echo "... scene.i"
	@echo "... scene.s"
	@echo "... shade.obj"
	@echo "... shade.i"
	@echo "... shade.s"
	@echo "... spectrum.obj"
	@echo "... spectrum.i"
	@echo "... spectrum.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

