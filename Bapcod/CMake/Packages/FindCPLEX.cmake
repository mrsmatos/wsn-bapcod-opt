# - Try to find Cplex 12.8, 12.9, 12.10, 20.0, 20.1, 22.1

if("${CPLEX_INCLUDE_DIR}" STREQUAL "" OR "${CPLEX_INCLUDE_DIR}" STREQUAL "CPLEX_INCLUDE_DIR-NOTFOUND")
	include(Solver)

	if(APPLE)
		set(Extension ".dylib")
	elseif(UNIX)
		set(Extension ".so")
	else()
		set(Extension ".lib")
	endif()

	if(EXISTS "$ENV{CPLEX_ROOT}/cplex/include/ilcplex")
		SET(CPLEX_INCLUDE_DIR "$ENV{CPLEX_ROOT}/cplex/include/ilcplex")
	endif()

	if(UNIX OR APPLE)
		if(EXISTS "$ENV{CPLEX_ROOT}/cplex/bin/")
			file(GLOB CPLEX_LIBRARY "$ENV{CPLEX_ROOT}/cplex/bin/*/libcplex*${Extension}")
		endif()
	else()
		if(EXISTS "$ENV{CPLEX_ROOT}/cplex/lib/")
			if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
				file(GLOB CPLEX_LIBRARY "$ENV{CPLEX_ROOT}/cplex/lib/*/stat_mdd/cplex*${Extension}")
			else()
				file(GLOB CPLEX_LIBRARY "$ENV{CPLEX_ROOT}/cplex/lib/*/stat_mda/cplex*${Extension}")
			endif()
		endif()
	endif ()

	include(FindPackageHandleStandardArgs)
	find_package_handle_standard_args(CPLEX DEFAULT_MSG CPLEX_INCLUDE_DIR CPLEX_LIBRARY)

	mark_as_advanced(CPLEX_LIBRARY CPLEX_INCLUDE_DIR)

endif()

if(NOT "${CPLEX_LIBRARY}" STREQUAL "CPLEX_LIBRARY-NOTFOUND")
	append_solver_libraries(${CPLEX_LIBRARY})
	append_solver_definitions("-D_CPLEX_FOUND")
endif()

if(NOT "${CPLEX_INCLUDE_DIR}" STREQUAL "CPLEX_INCLUDE_DIR-NOTFOUND")
	append_solver_include_dirs(${CPLEX_INCLUDE_DIR})
endif()