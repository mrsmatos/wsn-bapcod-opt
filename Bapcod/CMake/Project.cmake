# Macro to set a variable if it's not already set.
# Usage: set_if_not_set(VAR [value1] [value2]...)
macro(set_if_not_set name)
  if(NOT DEFINED "${name}")
      set(${name} ${ARGN})
  endif()
endmacro(set_if_not_set name value)

# Set the project name with the parent directory name in prefix.
# Usage: set_project_with_parent_name_prefix(a_project_name)
macro(set_project_with_parent_name_prefix PROJECT_NAME)
  #get_filename_component(parent_name ${CMAKE_CURRENT_SOURCE_DIR} PATH)
  #get_filename_component(parent_name ${parent_name} NAME)
  #string(TOLOWER ${parent_name} parent_name)
  project(${PROJECT_NAME} ${ARGN})
endmacro(set_project_with_parent_name_prefix PROJECT_NAME)

# Test if a variable is not empty. If the variable is empty then a FATAL_ERROR is throw.
# Usage: test_if_not_empty(VAR)
macro(test_if_not_empty name)
  if ("${${name}}" STREQUAL "")
    message(FATAL_ERROR "${name} is empty")
  endif()
endmacro(test_if_not_empty name)

# Get the name of the current source folder and set the project name with it.
# Usage: set_auto_project_name
macro(set_auto_project_name)
    get_filename_component(current_dir_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    set_project_with_parent_name_prefix(${current_dir_name} ${ARGN})
endmacro(set_auto_project_name)

# Initialize the project variables with default values if they are not already set.
# Usage: init_project_variables()
macro(init_project_variables)
    if ("${PROJECT_NAME}" STREQUAL "Project")
        set_auto_project_name()
    endif()

    # Variables to define application include/source directories, application output name and output directory.
    # Those variables can be overload if the user specify it inside the project's CMakeLists.txt.
    set_if_not_set(PROJECT_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
    set_if_not_set(PROJECT_SOURCE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
    set_if_not_set(PROJECT_OUTPUT_NAME "${PROJECT_NAME}")

    # Variables to other libraries
    # Those variables can be overload if the user specify it inside the project's CMakeLists.txt.
    set_if_not_set(ADDITIONAL_LIBRARIES "")
    set_if_not_set(ADDITIONAL_INCLUDE_DIRS "")
endmacro(init_project_variables )

# Search headers, source files and add include directories and link directories.
# Usage: search_and_add_files_to_build_the_project()
macro(search_and_add_files_to_build_the_project)
#    include_directories(${PROJECT_INCLUDE_DIRS} ${ADDITIONAL_INCLUDE_DIRS})
    include_directories(${ADDITIONAL_INCLUDE_DIRS})
    if(XCODE)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${Boost_INCLUDE_DIR}") #http://stackoverflow.com/questions/1867065/how-to-suppress-gcc-warnings-from-library-headers
      #message("BOOST directory ${Boost_INCLUDE_DIR} is also added as CMAKE_CXX_FLAGS -isystem") 
    else()
      include_directories(SYSTEM ${Boost_INCLUDE_DIR}) #This suppress the warnings rooted from Boost
      #message("BOOST directory ${Boost_INCLUDE_DIR} is also added as a SYSTEM") #http://stackoverflow.com/questions/1867065/how-to-suppress-gcc-warnings-from-library-headers
    endif()
    link_directories(${ADDITIONAL_LIBRARY_DIRS})

    if ("${header_files}" STREQUAL "")
      foreach(include_dir ${PROJECT_INCLUDE_DIRS})
	include_directories(${include_dir})
	file(GLOB tmp_headers
          ${include_dir}/*.hpp
          ${include_dir}/*.hxx
          ${include_dir}/*.h)
	list(APPEND header_files ${tmp_headers})
      endforeach()
    endif()

    if ("${source_files}" STREQUAL "")
      foreach(source_dir ${PROJECT_SOURCE_DIRS})
	file(GLOB tmp_sources
          ${source_dir}/*.cpp
          ${source_dir}/*.cxx
          ${source_dir}/*.c)
	list(APPEND source_files ${tmp_sources})
      endforeach()
    endif()

endmacro(search_and_add_files_to_build_the_project)

# Add target '${PROJECT_OUTPUT_NAME}' to build a static library.
# Usage: try_build_library()
macro(try_build_library)
    #message("Compiler id : ${CMAKE_CXX_COMPILER_ID} - CLANG_LIB_STD set to : ${LANG_LIB_STD}")

    init_project_variables()

    set_if_not_set(PROJECT_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/${CMAKE_BUILD_TYPE}")

    search_and_add_files_to_build_the_project()
    test_if_not_empty(source_files)
    add_library(${PROJECT_OUTPUT_NAME} STATIC ${source_files} ${header_files} ${ADDITIONAL_FILES})
    add_library(${PROJECT_OUTPUT_NAME}-shared SHARED ${source_files} ${header_files} ${ADDITIONAL_FILES})


    if(${Boost_VERSION} GREATER 107600 AND APPLE)
        find_package(ICU COMPONENTS  data  i18n iuuc )
        target_link_libraries(${PROJECT_OUTPUT_NAME} ${ICU_LIBRARIES})
        get_filename_component(ICU_LIBRARY_DIR ${ICU_LIBRARIES} DIRECTORY)
        target_link_directories(${PROJECT_OUTPUT_NAME} PUBLIC ${ICU_LIBRARY_DIR})
    endif()


    target_link_libraries(${PROJECT_OUTPUT_NAME} ${ADDITIONAL_LIBRARIES})

    target_link_libraries(${PROJECT_OUTPUT_NAME}-shared ${ADDITIONAL_LIBRARIES})

    if(UNIX)
      if("${CLP_FOUND}" STREQUAL "FOUND")
        target_link_libraries(${PROJECT_OUTPUT_NAME}-shared z readline bz2 ncurses blas lapack)
      endif()
    endif()

    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES PROJECT_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_LIBRARY_OUTPUT_DIRECTORY}"
                                                            ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_LIBRARY_OUTPUT_DIRECTORY}")

endmacro(try_build_library)

# Add target '${PROJECT_OUTPUT_NAME}' to build a static library.
# Usage: try_build_shared_library()
macro(try_build_shared_library)

    init_project_variables()

    set_if_not_set(PROJECT_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/${CMAKE_BUILD_TYPE}")

    search_and_add_files_to_build_the_project()
    test_if_not_empty(source_files)
    add_library(${PROJECT_OUTPUT_NAME} SHARED ${source_files} ${header_files} ${ADDITIONAL_FILES})

    target_link_libraries(${PROJECT_OUTPUT_NAME} ${ADDITIONAL_LIBRARIES})

    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES PROJECT_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_LIBRARY_OUTPUT_DIRECTORY}"
                                                            ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_LIBRARY_OUTPUT_DIRECTORY}")

endmacro(try_build_shared_library)

# Add target '${PROJECT_OUTPUT_NAME}' to build an application.
# Usage: try_build_application()
macro(try_build_application)

    init_project_variables()

    set_if_not_set(PROJECT_EXECUTABLE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE}")

    search_and_add_files_to_build_the_project()
    test_if_not_empty(source_files)
    add_executable(${PROJECT_OUTPUT_NAME} ${source_files} ${header_files} ${ADDITIONAL_FILES})

    target_link_libraries(${PROJECT_OUTPUT_NAME} ${ADDITIONAL_LIBRARIES})

    set_target_properties(${PROJECT_OUTPUT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_EXECUTABLE_OUTPUT_DIRECTORY})
    
endmacro(try_build_application)
