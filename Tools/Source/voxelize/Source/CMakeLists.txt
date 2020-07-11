
cmake_minimum_required(VERSION 2.8)

PROJECT(voxelize)

find_package(VTK 8.1 REQUIRED NO_MODULE)
#include(${VTK_USE_FILE})

add_executable(voxelize MACOSX_BUNDLE voxelize)

if(VTK_LIBRARIES)
  target_link_libraries(voxelize ${VTK_LIBRARIES})
else()
  target_link_libraries(voxelize vtkHybrid vtkWidgets)
endif()
