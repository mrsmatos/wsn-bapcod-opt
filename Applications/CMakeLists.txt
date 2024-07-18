cmake_minimum_required (VERSION 3.12)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/CMake)

add_subdirectories(
# put here your application(s)
BinAndVectorPacking
BinPackingAndCuttingStock
CapacitatedVehicleRouting
ColoredTSP
CVRP-Joao
GeneralizedAssignment
InventoryRouting
LotSizing
MultiCommodityFlow
Renault
VertexColoring
TimeDependentVehicleRoutingWithTimeWindows
VehicleRoutingWithDrones
)

# Try to make a package containing all demos.
set(PROJECT_OUTPUT_NAME "Applications")
set_project_version(3 0 0)
include(BapcodAllApplicationPackage)
try_build_bapcod_all_application_package()

macro(subdirlist result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/* ${curdir}/*/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
      string(REGEX REPLACE "\\//" "\\\\" child ${child})
      set(dirlist ${dirlist} ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

subdirlist(subdirs ${CMAKE_CURRENT_SOURCE_DIR})
