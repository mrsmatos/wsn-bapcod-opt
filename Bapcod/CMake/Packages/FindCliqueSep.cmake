################################################################
# This files permits to find CliqueSep source code
# (this code permits to use clique separation with BaPCod)
#
# * CLQSEP_INCLUDE_DIR: CliqueSep include dir
# * CLQSEP_SOURCE_DIR : CliqueSep source dir 
#
################################################################

if ("${CLQSEP_INCLUDE_DIR}" STREQUAL "" OR "${CLQSEP_INCLUDE_DIR}" STREQUAL "CLQSEP_INCLUDE_DIR-NOTFOUND")
  find_path(CLQSEP_INCLUDE_DIR
    NAMES clique.h
    PATHS ${CMAKE_SOURCE_DIR}/Tools/separation_libs/CliqueSep/include
    PATHS ENV{BAPCOD_ROOT}/Tools/separation_libs/CliqueSep/include) # for users apps

  find_path(CLQSEP_SOURCE_DIR 
    NAMES clique.cpp
    PATHS ${CMAKE_SOURCE_DIR}/Tools/separation_libs/CliqueSep/src
    PATHS ENV{BAPCOD_ROOT}/Tools/separation_libs/CliqueSep/src) # for users apps

  if (("${CLQSEP_INCLUDE_DIR}" STREQUAL "") OR ("${CLQSEP_INCLUDE_DIR}" STREQUAL "CLQSEP_INCLUDE_DIR-NOTFOUND"))
    	set(CLQSEP_INCLUDE_DIR "")
    	set(CLQSEP_SOURCE_DIR "")
  else()     
       add_definitions(-DCLQ_SEP_IS_FOUND)
       message("-- Found optional CliqueSep library")    
  endif()

  #include(FindPackageHandleStandardArgs)
  #find_package_handle_standard_args(CliqueSep "Could not find optional CliqueSep lib" CLQSEP_INCLUDE_DIR CLQSEP_SOURCE_DIR)
  mark_as_advanced(CLQSEP_SOURCE_DIR CLQSEP_INCLUDE_DIR)
endif()
