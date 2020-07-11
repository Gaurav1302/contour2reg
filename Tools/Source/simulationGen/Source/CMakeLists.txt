
cmake_minimum_required(VERSION 2.8)

PROJECT(simulationGen)

find_package(VTK 8.1 REQUIRED NO_MODULE)
#include(${VTK_USE_FILE})

add_executable(simulationGen MACOSX_BUNDLE simulationGen)

if(VTK_LIBRARIES)
  target_link_libraries(simulationGen ${VTK_LIBRARIES})
else()
  target_link_libraries(simulationGen vtkHybrid vtkWidgets)
endif()
