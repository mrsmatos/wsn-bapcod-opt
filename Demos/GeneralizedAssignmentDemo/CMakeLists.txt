cmake_minimum_required (VERSION 3.12) 

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Bapcod/CMake
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Bapcod/CMake/Packages
    ${CMAKE_CURRENT_SOURCE_DIR}/CMake ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Packages)

include(BapcodApplicationProject)
include(BapcodApplicationPackage)
include(Documentation)
include(Profile)

# Optional: Set other libraries
# set(USER_LIBRARIES other_libraries)
# set(USER_INCLUDE_DIR other_include_libraries_dir)

set_project_version(1 0 0)

# Build Application
try_build_bapcod_application()

# Make Package
try_build_bapcod_application_package()

