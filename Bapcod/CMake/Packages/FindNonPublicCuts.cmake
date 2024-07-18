# try to find installed version

find_path(NON_PUBLIC_CUTS_INCLUDE_DIR
	  NAMES bcModelNonPublicCuts.hpp
	  PATHS "${CMAKE_SOURCE_DIR}/Bapcod/include"
	  )

if (NOT "${NON_PUBLIC_CUTS_INCLUDE_DIR}" STREQUAL "NON_PUBLIC_CUTS_INCLUDE_DIR-NOTFOUND")
    add_definitions(-DUSE_NON_PUBLIC_CUTS)
    message("-- Found optional non-public cuts sources")    
endif()
