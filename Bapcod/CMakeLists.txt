cmake_minimum_required (VERSION 3.12)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/CMake $ENV{BAPCOD_ROOT}/CMake/Packages ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Packages ${CMAKE_CURRENT_SOURCE_DIR}/../Tools/lemon-1.3.1/lemon)
#message(CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH})

set(BAPCOD_VERSION 0.82.5 CACHE INTERNAL "")
set(BAPCOD_VERSION_DATE 24/09/2023 CACHE INTERNAL "")

#Project Name
project(bapcod VERSION ${BAPCOD_VERSION} LANGUAGES CXX)

set(CLANG_LIB_STD OFF CACHE BOOL "clang needs to use libstd")
if(CLANG_LIB_STD)
  message( "CLANG_LIB_STD set to : " ${CLANG_LIB_STD} )
endif()


if (WIN32)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(LEMON_LIBRARY "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/lemon-1.3.1/build/lemon/Debug/lemon.lib" CACHE FILEPATH "LEMON library" FORCE )
    endif()   
endif()

include(BapcodProject)
include(BapcodPackage)
include(Documentation)
include(Profile)
include(BapcodTest)

# Build Bapcod
try_build_bapcod()

# Make Package
try_build_bapcod_package()

# Build Documentation
# set(${PROJECT_OUTPUT_NAME}_DOC_EXCLUDE_PATTERNS "*/include_dev/* */include_utils/*")
#try_build_documentation()

# Make Profile
try_build_profile()

# Make unit test with cxxTest
try_build_bapcod_test()
