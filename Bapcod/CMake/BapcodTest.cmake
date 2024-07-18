include(Test)

macro(try_build_bapcod_test)

  set(BUILD_TESTS "OFF" 
    CACHE BOOL "Activate/Deactivate unit_test building")

  if (BUILD_TESTS)
    set(TEST_NAME ${PROJECT_OUTPUT_NAME}_programParametersTest)
    set(TEST_INPUT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test/unit/cxxtest/include)
    set(TEST_INPUT_INCLUDE_DIRS ${TEST_INPUT_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/test/common/include ${APPLICATION_INCLUDE_DIRS})
    set(TEST_OUTPUT_CPP_FILE ProgramParameters.cpp)

    file(GLOB_RECURSE TEST_INPUT_FILES
      ${TEST_INPUT_INCLUDE_DIR}/*.hpp
      ${TEST_INPUT_INCLUDE_DIR}/*.h
      ${TEST_INPUT_INCLUDE_DIR}/*.hxx)
  
    set(TEST_LIBRARY_INCLUDE_DIRS ${ADDITIONAL_LIBRARY_DIRS})

    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test/unit/cxxtest/config DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    try_build_test("${TEST_NAME}" "${TEST_INPUT_INCLUDE_DIRS}" "${TEST_OUTPUT_CPP_FILE}" "${TEST_INPUT_FILES}" "${TEST_LIBRARY_INCLUDE_DIRS}" "${TEST_LIBRARIES}")
  endif()

endmacro(try_build_bapcod_test)
