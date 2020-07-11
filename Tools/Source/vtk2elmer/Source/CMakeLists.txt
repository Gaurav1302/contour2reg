
cmake_minimum_required(VERSION 2.8)

PROJECT(vtk2elmer)

find_package(VTK 8.1 REQUIRED NO_MODULE)
#include(${VTK_USE_FILE})

add_executable(vtk2elmer MACOSX_BUNDLE vtk2elmer writeElmerFiles)

if(VTK_LIBRARIES)
  target_link_libraries(vtk2elmer ${VTK_LIBRARIES})
else()
  target_link_libraries(vtk2elmer vtkHybrid vtkWidgets)
endif()
