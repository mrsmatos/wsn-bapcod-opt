include(Solver)

if(NOT EXISTS $ENV{CLP_ROOT} )
  set(CLP_ROOT $ENV{COINOR_ROOT})
else()
  set(CLP_ROOT $ENV{CLP_ROOT})
endif ()

if("${CLP_INCLUDE_DIR}" STREQUAL "" OR "${CLP_INCLUDE_DIR}" STREQUAL "CLP_INCLUDE_DIR-NOTFOUND")

  SET (CLP_FOUND "NOT-FOUND" CACHE STRING "Found CLP Solver" FORCE)

  find_path(CLP_INCLUDE_DIR
    NAMES ClpSolve.hpp
    PATHS "${CLP_ROOT}/include/*"
    PATHS "${CLP_ROOT}/src/")

  find_library(CLP_LIBRARY
    NAMES libClp.lib libClp.a libClp.so libClp.dylib
    PATHS "${CLP_ROOT}/lib"
    PATHS "C:/Program Files/CoinCLP/MSVisualStudio/v16/x64/${CMAKE_BUILD_TYPE}")

  find_library(OSICLP_LIBRARY
    NAMES libOsiClp.lib libOsiClp.a libOsiClp.so libOsiClp.dylib
    PATHS "${CLP_ROOT}/lib"
    PATHS "C:/Program Files/CoinCLP/MSVisualStudio/v16/x64/${CMAKE_BUILD_TYPE}")

  find_path(COIN_INCLUDE_DIR
    NAMES CoinPragma.hpp
    PATHS "${CLP_ROOT}/include/*"
    PATHS "C:/Program Files/CoinUtils/src")

  find_library(COIN_LIBRARY
    NAMES libCoinUtils.lib libCoinUtils.a libCoinUtils.so libCoinUtils.dylib
    PATHS "${CLP_ROOT}/lib"
    PATHS "C:/Program Files/CoinCLP/MSVisualStudio/v16/x64/${CMAKE_BUILD_TYPE}")

 
 
  #include(FindPackageHandleStandardArgs)
  #find_package_handle_standard_args(CLP DEFAULT_MSG CLP_LIBRARY CLP_INCLUDE_DIR COIN_LIBRARY COIN_INCLUDE_DIR)

  if (NOT "${CLP_LIBRARY}" STREQUAL "CLP_LIBRARY-NOTFOUND" AND NOT "${COIN_LIBRARY}" STREQUAL "COIN_LIBRARY-NOTFOUND"
          AND  NOT "${COIN_INCLUDE_DIR}" STREQUAL "COIN_INCLUDE_DIR-NOTFOUND"
          AND  NOT "${CLP_INCLUDE_DIR}" STREQUAL "CLP_INCLUDE_DIR-NOTFOUND")
    SET (CLP_FOUND "FOUND" CACHE STRING "Found CLP Solver" FORCE)
    message("-- Found CLP solver")
  endif()

  mark_as_advanced(CLP_FOUND CLP_INCLUDE_DIR CLP_LIBRARY CLP_LIBRARY_DIR COIN_INCLUDE_DIR COIN_LIBRARY COIN_LIBRARY_DIR)

endif()


if (NOT "${CLP_LIBRARY}" STREQUAL "CLP_LIBRARY-NOTFOUND" AND NOT "${COIN_LIBRARY}" STREQUAL "COIN_LIBRARY-NOTFOUND"
    AND  NOT "${COIN_INCLUDE_DIR}" STREQUAL "COIN_INCLUDE_DIR-NOTFOUND"
    AND  NOT "${CLP_INCLUDE_DIR}" STREQUAL "CLP_INCLUDE_DIR-NOTFOUND")
  get_filename_component(CLP_LIBRARY_DIR "${CLP_LIBRARY}" PATH)
  get_filename_component(COIN_LIBRARY_DIR "${COIN_LIBRARY}" PATH)
  append_solver_include_dirs(${CLP_INCLUDE_DIR} ${COIN_INCLUDE_DIR})
  append_solver_libraries( ${OSICLP_LIBRARY} ${CLP_LIBRARY} ${COIN_LIBRARY})
  append_solver_library_dirs(${CLP_LIBRARY_DIR} ${COIN_LIBRARY_DIR})
  append_solver_definitions("-D_CLP_FOUND")

endif()
