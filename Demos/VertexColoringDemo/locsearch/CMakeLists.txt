# This file is not a part of bapcod project, it could not be developped in the bapcod's cmake environment

cmake_minimum_required (VERSION 3.12)
project (vcp_demo_locsearch C)
 
file (
	GLOB
	VCP_DEMO_LS_SOURCE_FILES
	*.c
	*.h
)

# The library will be generated in repertory lib
set(LIBRARY_OUTPUT_PATH "lib")

set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}")

# add the library
add_library (vcp_demo_locsearch ${VCP_DEMO_LS_SOURCE_FILES})

# Write the path to Sewell code in the cache
set(VCP_DEMO_LOCSEARCH_INCLUDE_DIR "${PROJECT_SOURCE_DIR}" CACHE INTERNAL "This is the path for include files of the local search algorithm by Held, Cook, and Sewell")
