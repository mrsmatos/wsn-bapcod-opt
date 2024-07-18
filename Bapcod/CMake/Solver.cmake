# - Set the Solver build options for all projects

include(Project)

macro(append_solver_libraries lib)
    set_if_not_set(SOLVER_LIBRARIES "")
    list(APPEND SOLVER_LIBRARIES ${lib} ${ARGN})
endmacro()

macro(append_solver_library_dirs lib_dir)
    set_if_not_set(SOLVER_LIBRARY_DIRS "")
    list(APPEND SOLVER_LIBRARY_DIRS ${lib_dir} ${ARGN})
endmacro()

macro(append_solver_include_dirs include_dir)
    set_if_not_set(SOLVER_INCLUDE_DIRS "")
    list(APPEND SOLVER_INCLUDE_DIRS ${include_dir} ${ARGN})
endmacro()

macro(append_solver_definitions def)
    set_if_not_set(SOLVER_DEFINITIONS "")
    list(APPEND SOLVER_DEFINITIONS ${def} ${ARGN})
endmacro()

set(SUPPORTED_SOLVERS "CPLEX;CLP")

# Search for all default solver's library and include dir
macro(try_find_solvers)
    foreach(SOLVER_VAR ${SUPPORTED_SOLVERS})
        # Upper and lower case name of the solver
        string(TOUPPER "${SOLVER_VAR}" SOLVER_UPPER)

        # Search for used solver's library and include directory
        find_package(${SOLVER_VAR})
    endforeach()
    
    # Upper case name of all the solvers, usefull if user redefine it
    string(TOUPPER "${SUPPORTED_SOLVERS}" UPPER_SUPPORTED_SOLVERS)
endmacro(try_find_solvers)
