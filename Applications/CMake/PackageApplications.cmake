# Set the version for the project
# Usage: set_project_name_version(MAJOR MINOR PATCH)
macro(set_project_version VERSION_MAJOR VERSION_MINOR VERSION_PATCH)
    string(TOUPPER "${PROJECT_NAME}" PROJECT_OUTPUT_NAME_UPPER_TMP)
    set_if_not_set(PROJECT_OUTPUT_NAME_UPPER "${PROJECT_OUTPUT_NAME_UPPER_TMP}")
    set_if_not_set(${PROJECT_OUTPUT_NAME_UPPER}_VERSION_MAJOR ${VERSION_MAJOR})
    set_if_not_set(${PROJECT_OUTPUT_NAME_UPPER}_VERSION_MINOR ${VERSION_MINOR})
    set_if_not_set(${PROJECT_OUTPUT_NAME_UPPER}_VERSION_PATCH ${VERSION_PATCH})
    set_if_not_set(${PROJECT_OUTPUT_NAME_UPPER}_VERSION "${${PROJECT_OUTPUT_NAME_UPPER}_VERSION_MAJOR}${${PROJECT_OUTPUT_NAME_UPPER}_VERSION_MINOR}${${PROJECT_OUTPUT_NAME_UPPER}_VERSION_PATCH}")
    
endmacro()


# Add target '${PROJECT_OUTPUT_NAME}_package' to build binary package of the project.
# Add target '${PROJECT_OUTPUT_NAME}_package_source' to build source package of the project.
# Add target 'package' to build all binary packages.
# Add target 'package_source' to build all source packages.
# Usage: try_build_package()
macro(try_build_package)
  set_project_version(1 0 0)
  
  set(CPACK_PACKAGE_NAME "${PROJECT_OUTPUT_NAME}")
  set(CPACK_PACKAGE_VENDOR "Romain LEGUAY")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Example test framework")
  set(CPACK_PACKAGE_VERSION_MAJOR "${${PROJECT_OUTPUT_NAME_UPPER}_VERSION_MAJOR}")
  set(CPACK_PACKAGE_VERSION_MINOR "${${PROJECT_OUTPUT_NAME_UPPER}_VERSION_MINOR}")
  set(CPACK_PACKAGE_VERSION_PATCH "${${PROJECT_OUTPUT_NAME_UPPER}_VERSION_PATCH}")
  set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
  set(CPACK_SOURCE_IGNORE_FILES "/lib/;/.svn/;/.git/;/.project/;/.cproject/;/CMakeFiles/;/*#/;/CMakeFiles/;/cmake_install.cmake/;/CMakeCache.txt/;/CPackConfig.cmake/;/CPackSourceConfig.cmake/;/Makefile/;/*.sln/;/.DS_Store/;/.gitignore/;")
  set(CPACK_SOURCE_INSTALLED_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR};.;${CMAKE_CURRENT_SOURCE_DIR}/include;./include;${CMAKE_CURRENT_SOURCE_DIR}/src;./src;${CMAKE_CURRENT_SOURCE_DIR}/CMake;./CMake") 

  set(CPACK_PACKAGE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/packages")
  
  if (NOT ("${CMAKE_BUILD_TYPE}" STREQUAL ""))
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_BUILD_TYPE}")
  endif()
  
  set(CPACK_INSTALL_CMAKE_PROJECTS)
  foreach(component ${COMPONENT_NAMES})
    list(APPEND CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_CURRENT_BINARY_DIR};${PROJECT_OUTPUT_NAME};${component};/")
  endforeach()
  
  set(CPACK_OUTPUT_CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/CPackConfig.cmake")
  set(CPACK_SOURCE_OUTPUT_CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/CPackSourceConfig.cmake")
  
  add_custom_target(${PROJECT_OUTPUT_NAME}_package "cpack" "--config" "${CPACK_OUTPUT_CONFIG_FILE}")
  add_custom_target(${PROJECT_OUTPUT_NAME}_package_source "cpack" "--config" "${CPACK_SOURCE_OUTPUT_CONFIG_FILE}")
  
  # If the target 'package' doesn't exist, we create it
#  if (NOT TARGET "package")
#    add_custom_target(package)
#  endif()
  
  # Add the dependency ${PROJECT_OUTPUT_NAME}_package to the target 'package'
#  add_dependencies(package ${PROJECT_OUTPUT_NAME}_package)
  
  # If the target 'package_source' doesn't exist, we create it
#  if (NOT TARGET "package_source")
#    add_custom_target(package_source)
#  endif()
  
  # Add the dependency ${PROJECT_OUTPUT_NAME}_package_source to the target 'package_source'
#  add_dependencies(package_source ${PROJECT_OUTPUT_NAME}_package_source)
  
  include(CPack)
  
endmacro(try_build_package)
