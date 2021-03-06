A how to guide to installing the required software to run the KP Solver that can be downloaded through the following steps (Cmake Build Tool 3.8, Doxygen and Boost 1.64)

Starting from scratch this is what you need to do (no CLION) ...

# Virtual Machine
!!! virtual box requires Intel Virtualization Technology in BIOS to be enabled
Power up computer
During boot enter BIOS setup (use of of those keys to get into BIOS: DEL, F1, F2, F8, F10, F12)
Set "Intel Virtualization Technology" Enabled

Firefox: https://www.virtualbox.org/wiki/Linux_Downloads
VirtualBox/virtualbox-5.1_5.1.2-108956~Ubuntu~xenial_amd64.deb

# Ubuntu
Firefox: http://www.ubuntu.com/desktop 
ubuntu-16.04.1-desktop-amd64.iso

$: sudo apt-get update
$: sudo apt-get upgrade

# C++ Compiler
$: sudo apt-get install build-essential 

# Cmake Build Tool 3.8
!!! boost version 1.64 requires cmake version 3.8, therefore building cmake by hand
delete previous installation of cmake if installed

$: sudo apt-get purge cmake
Firefox: https://cmake.org/files/
Firefox: Click v3.8
Firefox: Click “cmake-3.8.2.tar.gz”
Firefox: Save
$: cd ~
$: mdkir sw
$: cd sw
$: mdkir cmake
$: cd cmake
$: tar xzvf ~/Downloads/cmake-3.8.2.tar.gz
$: cd cmake-3.8.2
$: ./bootstrap
$: make 
$: sudo make install
$: cmake --version
cmake version 3.8.2
Take a look at the tutorial: https://cmake.org/cmake-tutorial/

# Doxygen
$: sudo apt-get install doxygen
$: sudo apt-get install doxygen-doc
$: sudo apt-get install doxygen-gui
$: sudo apt-get install doxygen-latex
$: sudo apt-get install graphviz

# Latex
$: sudo apt-get install texlive-full
$: sudo apt-get install texmaker

# Boost 1.64
Firefox: http://www.boost.org/users/history/ 
Firefox: Click Version 1.64.0 "Download"
Firefox: Click “boost_1_64_0.tar.bz2”
Firefox: Save
$: cd ~/sw
$: mdkir boost
$: cd boost
$: tar --bzip2 -xf ~/Downloads/boost_1_64_0.tar.bz2
$: cd boost_1_64_0
$: ./bootstrap.sh
$: sudo ./b2 install
# - Program Option Package
Take a look at: https://theboostcpplibraries.com/boost.program_options

"CMakeLists.txt" 

set(Boost_USE_STATIC_LIBS ON)
FIND_PACKAGE(Boost 1.64.0 COMPONENTS program_options filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(kp main.cpp)
target_link_libraries(kp ${Boost_LIBRARIES})

# Project
$: cd ~/sw/project

$: doxygen -g
$: gedit Doxyfile
PROJECT_NAME = "KP"
PROJECT_NUMBER = 1.0.0
EXTRACT_ALL = YES
EXTRACT_PRIVATE = YES
BUILTIN_STL_SUPPORT = YES
$: doxygen
# HTML
Double Click on ./html/index.html
# Latex
./latex

$: mkdir target
$: cd target 
$: cmake -DCMAKE_BUILD_TYPE=Release ..
# during debugging use: cmake -DCMAKE_BUILD_TYPE=Debug ..
$: sudo make install
/usr/local/bin/kp
# You can now use the progam kp in any directory
$: cd ~/data/thesis
$: KP --help
$: KP

# Create a deb package for proper distribution
$: cd ~/sw/project
$: mkdir release
$: cd release  
$: cmake ..
$: cpack
~/sw/project/release/kp-1.0.0-Linux.deb
# You can now distribute and install the program on any computer (Linux)
Copy ~/sw/project/release/kp-1.0.0-Linux.deb to the computer
$: sudo dpkg -i kp-1.0.0-Linux.deb
$: dpkg-query -l "kp"
# To uninstall the program
$: sudo apt-get remove kp
