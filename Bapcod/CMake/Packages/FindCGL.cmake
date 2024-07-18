
find_path(CGL_INCLUDE_DIR
  NAMES ClpSolve.hpp
  PATHS "$ENV{COINOR_ROOT}/include/coin"
  PATHS "/usr/local/include/coin"
  PATHS "/usr/include/coin"
  PATHS "C:/Program Files/CoinCLP/include")

find_library(CLP_LIBRARY
  NAMES Clp
  PATHS "$ENV{COINOR_ROOT}/lib"
  PATHS "/usr/local/lib/coin"
  PATHS "C:/Program Files/CoinCLP/lib")

find_library(CGL_LIBRARY
  NAMES Cgl
  PATHS "$ENV{COINOR_ROOT}/lib"
  PATHS "/usr/local/lib/coin"
  PATHS "C:/Program Files/CoinCLP/lib")

find_library(OSICLP_LIBRARY
  NAMES OsiClp
  PATHS "$ENV{COINOR_ROOT}/lib"
  PATHS "/usr/local/lib/coin"
  PATHS "C:/Program Files/CoinCLP/lib")

find_library(OSI_LIBRARY
  NAMES Osi
  PATHS "$ENV{COINOR_ROOT}/lib"
  PATHS "/usr/local/lib/coin"
  PATHS "C:/Program Files/CoinCLP/lib")

#    find_library(OSICBC_LIBRARY
#    NAMES OsiCbc
#    PATHS "/usr/local/lib/coin"
#    PATHS "C:/Program Files/CoinCLP/lib")

find_library(COIN_LIBRARY
  NAMES CoinUtils
  PATHS "$ENV{COINOR_ROOT}/lib"
  PATHS "/usr/local/lib/coin"
  PATHS "C:/Program Files/CoinCLP/lib")

#include(FindPackageHandleStandardArgs)
#find_package_handle_standard_args(CGL "Could not find optional CGL lib" CLP_LIBRARY CGL_INCLUDE_DIR COIN_LIBRARY)

mark_as_advanced(CGL_INCLUDE_DIR CGL_LIBRARY CLP_LIBRARY CLP_LIBRARY_DIR COIN_LIBRARY COIN_LIBRARY_DIR)

#message("CGL_LIBRARY     : ${CGL_LIBRARY}")
#message("CGL_INCLUDE_DIR : ${CGL_INCLUDE_DIR}")

set_if_not_set(CGL_INCLUDE_DIRS "")
set_if_not_set(CGL_LIBRARIES "")
set_if_not_set(CGL_LIBRARY_DIRS "")
set_if_not_set(CGL_DEFINITIONS "")

if (NOT (("${CLP_LIBRARY}" STREQUAL "CLP_LIBRARY-NOTFOUND") OR ("${COIN_LIBRARY}" STREQUAL "COIN_LIBRARY-NOTFOUND") OR ("${CGL_INCLUDE_DIR}" STREQUAL "CGL_INCLUDE_DIR-NOTFOUND")))
  get_filename_component(CLP_LIBRARY_DIR "${CLP_LIBRARY}" PATH)
  get_filename_component(COIN_LIBRARY_DIR "${COIN_LIBRARY}" PATH)
  
  list(APPEND CGL_INCLUDE_DIRS ${CGL_INCLUDE_DIR} )
  list(APPEND CGL_LIBRARIES ${OSI_LIBRARY} ${CGL_LIBRARY} ${OSICLP_LIBRARY} ${CLP_LIBRARY} ${COIN_LIBRARY} )
  list(APPEND CGL_LIBRARY_DIRS ${CLP_LIBRARY_DIR} ${COIN_LIBRARY_DIR} )
  list(APPEND CGL_DEFINITIONS "-D_CGL_FOUND" )

  message("-- Found optional CGL library")    
endif()

