cmake_minimum_required(VERSION 2.6)
project(raytr)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -march=native -DHAVE_SSE2 -DDSFMT_MEXP=19937")

include( PrecompiledHeader )

include_directories(SYSTEM "ext")

if(NOT NO_CV)
find_package( OpenCV REQUIRED core highgui imgproc imgcodecs )
else()
add_definitions(-DNO_CV)
endif()

set(SOURCE_FILES ext/json.hpp main.cpp camera.h common.h scene.cpp scene.h objekt.h primitives.h
    image.h primitives.cpp pathtracer.cpp pathtracer.h random.cpp random.h
    image.cpp photon_tree.cpp photon_tree.h pm.cpp pm.h shade.h ext/stb_image_write.h
    collisionmanager.h bbox_tree.cpp bbox_tree.h ext/SimpleObject.cpp ext/SimpleObject.h
    stdafx.h ext/tiny_obj_loader.cc ext/tiny_obj_loader.h objloader.h objloader.cpp
    shade.cpp bxdfs.h objekt.cpp spectrum.cpp spectrum.h
    ext/dSFMT/dSFMT.h ext/dSFMT/dSFMT.cpp meshsimp1.cpp meshsimp1.h)

add_executable(raytr ${SOURCE_FILES})
add_precompiled_header(raytr stdafx.h)
target_link_libraries(raytr ${OpenCV_LIBRARIES})
