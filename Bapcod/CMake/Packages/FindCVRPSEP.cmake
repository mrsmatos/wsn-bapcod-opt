#message("CVRPSEP searching in  : " ${CMAKE_CURRENT_SOURCE_DIR}/../Tools/separation_libs/CVRPSEP )

  find_path(CVRPSEP_INCLUDE_DIR
    NAMES basegrph.h
    PATHS ${CMAKE_SOURCE_DIR}/Tools/separation_libs/CVRPSEP/include # for dev apps
    PATHS ENV{BAPCOD_ROOT}/Tools/separation_libs/CVRPSEP/include) # for users apps

  find_path(CVRPSEP_SOURCE_DIR
    NAMES basegrph.cpp
    PATHS ${CMAKE_SOURCE_DIR}/Tools/separation_libs/CVRPSEP/src # for dev apps
    PATHS $ENV{BAPCOD_ROOT}/Tools/separation_libs/CVRPSEP/src) # for users apps

  if (("${CVRPSEP_INCLUDE_DIR}" STREQUAL "") OR ("${CVRPSEP_INCLUDE_DIR}" STREQUAL "CVRPSEP_INCLUDE_DIR-NOTFOUND"))
    	set(CVRPSEP_INCLUDE_DIR "")
    	set(CVRPSEP_SOURCE_DIR "")
  else()     
       add_definitions(-DCVRPSEP_IS_FOUND)
       message("-- Found optional CVRPSEP library")    
  endif()

  #include(FindPackageHandleStandardArgs)
  #find_package_handle_standard_args(CVRPSEP "Could not find optional CVRPSEP lib" CVRPSEP_SOURCE_DIR CVRPSEP_INCLUDE_DIR)

  mark_as_advanced(CVRPSEP_SOURCE_DIR CVRPSEP_INCLUDE_DIR)
