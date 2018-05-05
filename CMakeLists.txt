cmake_minimum_required(VERSION 3.8)
project(KP)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS_RELEASE "-s -O3 -mavx -msse2 -mfma -mbmi -Wno-unused-variable")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g")

# boost C++ libraries

set(Boost_USE_STATIC_LIBS ON)
FIND_PACKAGE(Boost 1.64.0 COMPONENTS program_options filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# program
add_executable(KP main.cpp Soliton.cpp Lump.cpp KP.cpp csv.cpp)
target_link_libraries(KP ${Boost_LIBRARIES})

install(TARGETS KP DESTINATION bin)

# deployment
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Athina Lange") # required
include(CPack)
