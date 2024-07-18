include(Project)

macro(copy_data_and_config_folders_to_build_folder)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/config DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endmacro(copy_data_and_config_folders_to_build_folder)

macro(create_symlinks_to_benchmarks_config_data_tests_folders)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/ ${CMAKE_CURRENT_BINARY_DIR}/benchmarks)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/config/ ${CMAKE_CURRENT_BINARY_DIR}/config)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data/ ${CMAKE_CURRENT_BINARY_DIR}/data)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/tests/ ${CMAKE_CURRENT_BINARY_DIR}/tests)
endmacro(create_symlinks_to_benchmarks_config_data_tests_folders)

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

macro(try_build_bapcod_application)

  set_auto_project_name(CXX)

  IF ( NOT WIN32 )
    use_cxx11()
  ENDIF()

  set(Boost_USE_STATIC_LIBS ON)
  find_package(Boost COMPONENTS program_options regex system timer chrono thread REQUIRED)
	
  find_package(Bapcod)
  #find_package(Backtrace)
  find_package(LEMON)
  find_package(RapidJSON)
  find_package(MathProgSolverBuilder)
  include(Solver)
    
  if("${SOLVER_LIBRARIES}" STREQUAL "")
    try_find_solvers()
  endif()

  # Set all the include dirs used by Bapcod and the current Application
  set(PROJECT_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${LEMON_INCLUDE_DIR}
    ${BAPCOD_INCLUDE_DIRS} 
    ${SOLVER_INCLUDE_DIRS} 
    ${Boost_INCLUDE_DIR} 
    ${MPSB_INCLUDE_DIR}
    ${BACKTRACE_INCLUDE_DIR}
    ${USER_INCLUDE_DIR}

    )

  #set additional include rapidjson
  if ( NOT "${RAPIDJSON_INCLUDE_DIR}" STREQUAL "RAPIDJSON_INCLUDE_DIR-NOTFOUND")
    set(PROJECT_INCLUDE_DIRS  ${PROJECT_INCLUDE_DIRS} ${RAPIDJSON_INCLUDE_DIR})
  endif()

  # Set all the source dirs used by Bapcod and the current Application
  set(PROJECT_SOURCE_DIRS src ${MPSB_SOURCE_DIR})
  
  if( NOT CMAKE_BUILD_TYPE )
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
  endif()
  
  # if the target bapcod is found
  if(TARGET bapcod)
    # set additional libraries (boost and solvers) and build bapcod
    set(ADDITIONAL_LIBRARIES bapcod ${Boost_LIBRARIES} ${SOLVER_LIBRARIES} ${LEMON_LIBRARY} )
  endif()

  # set additional library dirs.
  set(ADDITIONAL_LIBRARY_DIRS ${Boost_INCLUDE_DIR} ${SOLVER_LIBRARY_DIRS} ) 

  # Add Definitions flags for the current project
  if(${ADLIB})
    add_definitions("-DADLIB")
  endif()

  add_definitions(${SOLVER_DEFINITIONS} ${RAPIDJSON_DEFINITIONS})


  if (WIN32)
    copy_data_and_config_folders_to_build_folder()
  else()
    create_symlinks_to_benchmarks_config_data_tests_folders()
  endif()

  set(PROJECT_EXECUTABLE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

  try_build_application()
  
  if (NOT TARGET bapcod)
    # set additional libraries (boost and solvers)
    target_link_libraries(${PROJECT_NAME} ${BAPCOD_LIBRARY} ${Boost_LIBRARIES} ${SOLVER_LIBRARIES} ${LEMON_LIBRARY}) 
  endif()
  
  # Define some specific defintion for macOS
  if(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "-framework IOKit -framework CoreFoundation -Wl,-no_compact_unwind" 
      CACHE STRING " Flags used by the linker." FORCE)
  endif()

  # Add those library for Boost::backtrace
  if(MSVC)
    target_link_libraries(${PROJECT_NAME} dbghelp)
  else()
    if(CMAKE_COMPILER_IS_GNUCXX)
      target_link_libraries(${PROJECT_NAME} dl)
     endif()
  endif()

  # Add pthread and z library for unix.
  if(UNIX)
      target_link_libraries(${PROJECT_NAME} pthread z)
      if("${CLP_FOUND}" STREQUAL "FOUND")
        target_link_libraries(${PROJECT_NAME} readline bz2 ncurses lapack blas)
      endif()
  endif()
      
  # Add user libraries (OPTIONAL)
  target_link_libraries(${PROJECT_NAME} ${USER_LIBRARIES})
        
endmacro(try_build_bapcod_application)
