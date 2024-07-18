
include(Project)

macro(use_cxx11)
  include(CheckCXXCompilerFlag)

  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=gnu++11" COMPILER_SUPPORTS_GNUXX11)

  if(COMPILER_SUPPORTS_GNUXX11 OR COMPILER_SUPPORTS_CXX11)
    #message("---- Setting CXX flags")
    if (CMAKE_VERSION VERSION_LESS "3.1")
      message("---- CMAKE_VERSION is less than 3.1. CXX flags will be set manually")
      if (COMPILER_SUPPORTS_GNUXX11)
        set (CMAKE_CXX_FLAGS "-std=gnu++11 ${CMAKE_CXX_FLAGS}")
      else() #COMPILER_SUPPORTS_CXX11
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
      endif ()
    else ()
      #message("---- CMAKE_CXX_STANDARD is set to 11 - Success")
      #message("Warning: If the operation system is Windows, this might not work as expected. Cmake code needs to be updated to handle this case because CMAKE_CXX_STANDARD always sets -std=gnu++11.")
      set (CMAKE_CXX_STANDARD 11)
    endif ()

  else()
      message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif()
      #message("")
endmacro(use_cxx11)



macro(try_build_bapcod)

  IF ( NOT WIN32)
    use_cxx11()
  ENDIF()  

  set(Boost_USE_STATIC_LIBS ON)

  if (WIN32)
      set(BOOST_LIBRARYDIR ${BOOST_ROOT}/stage/lib)
  endif()

  find_package(Boost COMPONENTS program_options regex system timer chrono thread REQUIRED)

  configure_file(${PROJECT_SOURCE_DIR}/include_dev/bcBapcodVersion.hpp.in ${PROJECT_SOURCE_DIR}/include_dev/bcBapcodVersion.hpp)
  #configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bcBapcodVersion.h.in ${CMAKE_CURRENT_BINARY_DIR}/bcBapcodVersion.h)

  find_package(BcpRcsp)
  find_package(Bapcod)
  find_package(NonPublicCuts)
  #find_package(Backtrace)
  find_package(MathProgSolverBuilder)
  find_package(CliqueSep)
  find_package(RheccSep)
  find_package(LEMON)
  #find_package(CGL)
  find_package(CVRPSEP)
  find_package(RapidJSON)
  # find_package(Libxml2)
  include(Solver)

  if("${SOLVER_LIBRARIES}" STREQUAL "")
    try_find_solvers()
  endif()


  set(PROJECT_INCLUDE_DIRS ${BAPCOD_INCLUDE_DIRS} ${BACKTRACE_INCLUDE_DIR} ${SOLVER_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}
    ${Boost_INCLUDE_DIR} ${MPSB_INCLUDE_DIR} ${LEMON_INCLUDE_DIR} ${CGL_INCLUDE_DIRS} ${BCP_RCSP_INCLUDE_DIR}
    ${CVRPSEP_INCLUDE_DIR} ${RHECCSEP_INCLUDE_DIR} ${CLQSEP_INCLUDE_DIR} )

#additional include dir rapidjson
if ( NOT "${RAPIDJSON_INCLUDE_DIR}" STREQUAL "RAPIDJSON_INCLUDE_DIR-NOTFOUND")
    set(PROJECT_INCLUDE_DIRS  ${PROJECT_INCLUDE_DIRS} ${RAPIDJSON_INCLUDE_DIR})
endif()
  # Set all the source dirs used by Bapcod

  set(PROJECT_SOURCE_DIRS ${MPSB_SOURCE_DIR} ${BAPCOD_SOURCE_DIR} ${BAPCOD_SOURCE_SOLVERS_DIR} ${BACKTRACE_SOURCE_DIR} ${CLQSEP_SOURCE_DIR} ${RHECCSEP_SOURCE_DIR}
      ${CVRPSEP_SOURCE_DIR})

  if("${BAPCOD_LIBRARY_RELEASE}" STREQUAL "BAPCOD_LIBRARY_RELEASE-NOTFOUND")
    #message("BAPCOD_LIBRARY_RELEASE WAS NOT FOUND")
    set(ADDITIONAL_LIBRARIES ${Boost_LIBRARIES} ${SOLVER_LIBRARIES} ${CGL_LIBRARIES} ${BCP_RCSP_LIBRARY} ${LEMON_LIBRARY})
  else()
    #message("BAPCOD_LIBRARY_RELEASE WAS FOUND")
    set(ADDITIONAL_LIBRARIES ${Boost_LIBRARIES} ${SOLVER_LIBRARIES} ${CGL_LIBRARIES} ${BCP_RCSP_LIBRARY} ${BAPCOD_LIBRARY_RELEASE} ${LEMON_LIBRARY})
  endif()

  # set additional library dirs.
  set(ADDITIONAL_LIBRARY_DIRS ${Boost_INCLUDE_DIR} ${SOLVER_LIBRARY_DIRS} ${CGL_LIBRARY_DIRS} ${BCP_RCSP_LIBRARY_DIRS})

  # Export the linked libraries
  set(EXPORTED_LIBRARIES_FILE "${CMAKE_BINARY_DIR}/dependent_libraries.cmake")
  file(WRITE ${EXPORTED_LIBRARIES_FILE} "Set(DEPENDENT_LIBRARIES \"${Boost_LIBRARIES};${LEMON_LIBRARY}\")")

  if( NOT CMAKE_BUILD_TYPE )
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
  endif()

  # if Bapcod source are found
  if(NOT BAPCOD_SOURCE_DIR-NOT_FOUND)
    set(PROJECT_OUTPUT_NAME bapcod)

    # On Windows or Unix sytem and not in Mac OS X

    if (NOT WIN32) 
      set(CMAKE_CXX_FLAGS "-fPIC ${CMAKE_CXX_FLAGS}")
    endif()

    if (WIN32 OR UNIX AND NOT APPLE)
      # If 64 bits mode
      if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    	set(PROJECT_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/${CMAKE_SYSTEM_NAME}-x64)
      else()
	    set(PROJECT_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/${CMAKE_SYSTEM_NAME}-x86)
      endif()
    else() # Apple System
        set(PROJECT_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/${CMAKE_SYSTEM_NAME})
    endif()
    if (WIN32)
        if(MSVC11)
            message("Contenu de CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
            set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /MP /d2Zi+ /Ob2 /Ot /Oi")
            set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MP")
            message("Contenu de CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
        endif()
        if(MSVC14)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D \"DEFINED_TO_STRING\" /D \"_SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS\" /MP /bigobj")
        endif()
    endif()

    if (NOT WIN32)
      add_compile_options("-std=c++0x")
    endif()

    try_build_library()


    set(TIMER "ON"
      CACHE BOOL "Activate/Desactivate timers used by Bapcod")

    # Activate/Desactivate Timer Option
    if(TIMER)
      set(TIMER_DEFINITION "-D_TIMER")
    else(TIMER)
      set(TIMER_DEFINITION "")
    endif()

    # Add Definitions flags for the current project
    add_definitions(${DEFINITIONS} ${TIMER_DEFINITION} ${SOLVER_DEFINITIONS} ${CGL_DEFINITIONS} ${RAPIDJSON_DEFINITIONS} )
    if(${ADLIB})
      add_definitions("-DADLIB")
    endif()

    # Define some specific defintion for macOS
    if(APPLE)
      set(CMAKE_EXE_LINKER_FLAGS "-framework IOKit -framework CoreFoundation -Wl,-no_compact_unwind"
        CACHE STRING " Flags used by the linker." FORCE)
    endif()

    # Set the library name prefix.
    # On Windows sytem
    if (WIN32)
        # Test if we use Visual Studio 2010
        if(MSVC10)
            set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES DEBUG_POSTFIX "_vc100_debug")
            set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES RELEASE_POSTFIX "_vc100_release")
        # Or test if we use Visual Studio 2012
        elseif(MSVC11)
			set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} " -v")
            set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES DEBUG_POSTFIX "_vc110_debug")
            set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES RELEASE_POSTFIX "_vc110_release")
         endif()
    else() # On Unix system and Mac OS X
       set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES DEBUG_POSTFIX "_debug")
       set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES RELEASE_POSTFIX "_release")
    endif()

    # Set the output library for all the build system.
    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_LIBRARY_OUTPUT_DIRECTORY})
    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PROJECT_LIBRARY_OUTPUT_DIRECTORY})
    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PROJECT_LIBRARY_OUTPUT_DIRECTORY})
    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL ${PROJECT_LIBRARY_OUTPUT_DIRECTORY})
    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${PROJECT_LIBRARY_OUTPUT_DIRECTORY})

    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_LIBRARY_OUTPUT_DIRECTORY})
    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${PROJECT_LIBRARY_OUTPUT_DIRECTORY})
    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${PROJECT_LIBRARY_OUTPUT_DIRECTORY})
    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_MINSIZEREL ${PROJECT_LIBRARY_OUTPUT_DIRECTORY})
    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO ${PROJECT_LIBRARY_OUTPUT_DIRECTORY})

  endif()

endmacro(try_build_bapcod)
