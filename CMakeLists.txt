cmake_minimum_required(VERSION 3.7.2)
project(s2d CXX)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(EIGEN3_INCLUDE_DIR ~/usr/eigen-3.3.9)

# Eigen3
#find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})



# Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})


# OpenMP
find_package(OpenMP)
if (OpenMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	link_libraries(${OpenMP_CXX_LIBRARIES})
endif(OpenMP_FOUND)

# My optimization flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
if (COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "-march=native -DNDEBUG")
else(COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
	set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG")
endif(COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
CHECK_CXX_COMPILER_FLAG("-Ofast" COMPILER_OPT_OFAST_SUPPORTED)
if (COMPILER_OPT_OFAST_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")
else(COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif(COMPILER_OPT_OFAST_SUPPORTED)



add_executable(s2d mesh.cpp pde.cpp export.cpp  main.cpp)




