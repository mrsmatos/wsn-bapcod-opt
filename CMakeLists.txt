cmake_minimum_required (VERSION 3.12)

project(BapcodFramework)

# Substitua essas variáveis pelos caminhos corretos se não estiverem sendo passadas via linha de comando
set(CPLEX_INCLUDE_DIR "/opt/ibm/ILOG/CPLEX_Studio2211/cplex/include")
include_directories(${CPLEX_INCLUDE_DIR})

set(CPLEX_LIBRARY "/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic/libcplex.a")
link_directories(${CPLEX_LIBRARY})

#set(CPLEX_LIBRARY_DIR "/opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic")
#link_directories(${CPLEX_LIBRARY_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -mfpmath=sse -march=native -flto -g -m64 -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++17 -flto=auto")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -ftree-vectorize -mfpmath=sse -march=native -flto -g -m64 -fPIC -fexceptions -DNDEBUG")


cmake_policy(SET CMP0111 NEW)

# Substitua 'main.cpp' e 'bapcod' pelos seus arquivos fonte e nome do alvo, respectivamente
#add_executable(mybapcod main.cpp)
#target_link_libraries(mybapcod ${CPLEX_LIBRARY})

set(BCP_RCSP_MIN_REQUIRED_VERSION 0.6.7)
option(ADLIB "Bapcod depends on ADLIB" OFF)
if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "arm*")
  set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE INTERNAL "" FORCE)
endif()


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_MACOSX_RPATH 1)
endif()

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

set(Boost_NO_WARN_NEW_VERSIONS ON)

if(DEFINED ENV{BOOST_ROOT})
  set(BOOST_ROOT $ENV{BOOST_ROOT})
  set(Boost_ROOT $ENV{BOOST_ROOT})
  #message("BOOST_ROOT and Boost_ROOT hint set to " ${BOOST_ROOT} " from BOOST_ROOT env variable ")
elseif ((EXISTS "$ENV{BAPCOD_ROOT}") OR (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Tools/boost_1_76_0"))
  set(Boost_NO_BOOST_CMAKE  ON)
  set(Boost_NO_SYSTEM_PATHS ON)  
  if(EXISTS "$ENV{BAPCOD_ROOT}/BaPCod/Tools/boost_1_76_0/build")
    set(BOOST_ROOT $ENV{BAPCOD_ROOT}/BaPCod/Tools/boost_1_76_0/build)
  elseif(EXISTS "$ENV{BAPCOD_ROOT}/Tools/boost_1_76_0/build")
    set(BOOST_ROOT $ENV{BAPCOD_ROOT}/Tools/boost_1_76_0/build)
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Tools/boost_1_76_0")
    if (WIN32)
      set(BOOST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/Tools/boost_1_76_0/)
    else()
      set(BOOST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/Tools/boost_1_76_0/build/)
    endif()
  endif()
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)


if (NOT WIN32)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -fvisibility-inlines-hidden ")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/CMake)

# This function add a list of subdirectories if their exists
function(add_subdirectories)
foreach(DIR ${ARGV})
  if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${DIR})
    add_subdirectory(${DIR})
  endif()
endforeach()
endfunction()

add_subdirectories(Tools
		   Bapcod
		   Demos
		   Applications
		   Tests
		  )

# Copy Scripts folder in the binary directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Scripts DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
