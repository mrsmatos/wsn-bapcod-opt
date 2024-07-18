################################################################
# This files permits to find RheccSep source code
# (this code permits to use rounded homogeneous extended capacity cut separation with BaPCod)
#
# * RHECCSEP_INCLUDE_DIR: RheccSep include dir
# * RHECCSEP_SOURCE_DIR : RheccSep source dir 
#
################################################################

if ("${RHECCSEP_INCLUDE_DIR}" STREQUAL "" OR "${RHECCSEP_INCLUDE_DIR}" STREQUAL "RHECCSEP_INCLUDE_DIR-NOTFOUND")
  find_path(RHECCSEP_INCLUDE_DIR
    NAMES CutGenerator.h
    PATHS ${CMAKE_SOURCE_DIR}/Tools/separation_libs/RheccSep/include # for dev apps
    PATHS ENV{BAPCOD_ROOT}/Tools/separation_libs/RheccSep/include) # for users apps

  find_path(RHECCSEP_SOURCE_DIR 
    NAMES CutGenerator.cpp
    PATHS ${CMAKE_SOURCE_DIR}/Tools/separation_libs/RheccSep/src # for dev apps
    PATHS ENV{BAPCOD_ROOT}/Tools/separation_libs/RheccSep/src) # for users apps

  if (("${RHECCSEP_INCLUDE_DIR}" STREQUAL "") OR ("${RHECCSEP_INCLUDE_DIR}" STREQUAL "RHECCSEP_INCLUDE_DIR-NOTFOUND"))
    	set(RHECCSEP_INCLUDE_DIR "")
    	set(RHECCSEP_SOURCE_DIR "")
  else()     
       add_definitions(-DRHECC_SEP_IS_FOUND)
       message("-- Found optional RheccSep library")    
  endif()

  #include(FindPackageHandleStandardArgs)
  #find_package_handle_standard_args(RheccSep "Could not find optional RheccSep lib" RHECCSEP_INCLUDE_DIR RHECCSEP_SOURCE_DIR)
  mark_as_advanced(RHECCSEP_SOURCE_DIR RHECCSEP_INCLUDE_DIR)
endif()
