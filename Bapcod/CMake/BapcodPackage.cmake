include(Package)

macro(try_build_bapcod_package)

  set(BAPCOD_ROOT_DIR "BaPCod")
  set(PACKAGE_SOURCE_IGNORE_FILES "/test/")
  
  # Install the library in BaPCod/lib
  install(TARGETS ${PROJECT_OUTPUT_NAME}
    LIBRARY
    ARCHIVE DESTINATION "${BAPCOD_ROOT_DIR}/lib"
    COMPONENT lib)

  # Install the tools (Backtrace lib...) in BaPCod/
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools
    DESTINATION "${BAPCOD_ROOT_DIR}/"
    COMPONENT tools)

  # Install the license, authors, and readme files in BaPCod/
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/bcBapcodVersion.h 
    ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.pdf
    ${CMAKE_CURRENT_SOURCE_DIR}/AUTHORS.txt
    ${CMAKE_CURRENT_SOURCE_DIR}/README.txt 
    DESTINATION "${BAPCOD_ROOT_DIR}/"
    COMPONENT library_header_version)

  # Install the minimum header (include and include_solvers)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include 
    DESTINATION "${BAPCOD_ROOT_DIR}/"
    COMPONENT library_headers)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include_solverinterfaces
    DESTINATION "${BAPCOD_ROOT_DIR}/"
    COMPONENT library_headers)
  
  # Set all the components to be installed and packed
  set(COMPONENTS_NAMES lib tools library_header_version library_headers)

  # Add two options to test if the package will contain dev include or util include.
  set(PACKAGE_DEV_INCLUDE NO CACHE BOOL "Add dev include inside the package?")
  set(PACKAGE_UTIL_INCLUDE NO CACHE BOOL "Add util include inside the package?")

  if(PACKAGE_DEV_INCLUDE)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include_dev
      DESTINATION "${BAPCOD_ROOT_DIR}/"
      COMPONENT library_headers)
    
  endif()

  if(PACKAGE_DEV_INCLUDE OR PACKAGE_UTIL_INCLUDE)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include_utils
      DESTINATION "${BAPCOD_ROOT_DIR}/"
      COMPONENT library_headers)
  endif()

  set(CPACK_PACKAGE_NAME "${BAPCOD_ROOT_DIR}")
  set(CPACK_PACKAGE_VENDOR "François VANDERBECK")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "BaPCod is a prototype code that solves Mixed Integer Programs (MIP) by application of a Dantzig-Wolfe reformulation technique. The reformulated problem is solved using a branch-and-price (column generation) algorithm.")
  
  # All those lines set the name of the package
  if (PACKAGE_DEV_INCLUDE AND PACKAGE_UTIL_INCLUDE)
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-Full-${CMAKE_SYSTEM_NAME}")
    if (NOT ("${CMAKE_BUILD_TYPE}" STREQUAL ""))
      set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-Full-${CMAKE_SYSTEM_NAME}-${CMAKE_BUILD_TYPE}")
    endif()
  elseif(PACKAGE_DEV_INCLUDE AND NOT PACKAGE_UTIL_INCLUDE)
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-Dev-${CMAKE_SYSTEM_NAME}")
    if (NOT ("${CMAKE_BUILD_TYPE}" STREQUAL ""))
      set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-Dev-${CMAKE_SYSTEM_NAME}-${CMAKE_BUILD_TYPE}")
    endif()
  elseif(NOT PACKAGE_DEV_INCLUDE AND PACKAGE_UTIL_INCLUDE)
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-Utils-${CMAKE_SYSTEM_NAME}")
    if (NOT ("${CMAKE_BUILD_TYPE}" STREQUAL ""))
      set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-Utils-${CMAKE_SYSTEM_NAME}-${CMAKE_BUILD_TYPE}")
    endif()
  endif()

  if (APPLE)
    set(CPACK_GENERATOR "DragNDrop")
  endif()

  set(CPACK_PACKAGE_DEFAULT_LOCATION "${CPACK_PACKAGE_NAME}")
  set(CPACK_PACKAGE_CONTACT "François Vanderbeck: fv@math.u-bordeaux1.fr")

  if(UNIX AND NOT APPLE)
    set(CPACK_PACKAGING_INSTALL_PREFIX "/opt")
    set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "BaPCod: a generic Branch-And-Price Code ${CPACK_PACKAGE_DESCRIPTION_SUMMARY}")
    set(CPACK_GENERATOR "DEB;STGZ;TGZ;RPM")
  endif()

  if(WIN32)
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_DEFAULT_LOCATION}")
    set(CPACK_NSIS_INSTALL_ROOT "D:")
    set(CPACK_NSIS_PACKAGE_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY}")
    set(CPACK_NSIS_MODIFY_PATH ON)
    set(BAPCOD_INCLUDE_PATH "D:\\${CPACK_NSIS_PACKAGE_NAME}\\${CPACK_PACKAGE_NAME}\\include" "D:\\${CPACK_NSIS_PACKAGE_NAME}\\${CPACK_PACKAGE_NAME}\\lib")
    #set (CPACK_NSIS_EXTRA_INSTALL_COMMANDS "Push \\\"PATH\\\"\nPush \\\"A\\\"\nPush \\\"HKLM\\\"\nPush \\\"D:\\\\BaPCod\\\\BaPCod\\\\include;D:\\\\BaPCod\\\\BaPCod\\\\lib\\\"\nCall EnvVarUpdate\nPop \$1")
  endif()

  set(CPACK_SOURCE_GENERATOR "TGZ" "ZIP")

  try_build_package()

endmacro(try_build_bapcod_package)
