cmake_minimum_required(VERSION 2.8)
project(seg)
if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

FIND_PACKAGE ( ITK REQUIRED)
INCLUDE( ${ITK_USE_FILE} )
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()


add_executable(AnisotropicMain AnisotropicMain.cxx )
target_link_libraries(AnisotropicMain ${ITK_LIBRARIES})

add_executable(ADLBRMain ADLBRMain.cxx )
target_link_libraries(ADLBRMain ${ITK_LIBRARIES})


add_executable(OtherDenoisingMain OtherDenoisingMain.cxx )
target_link_libraries(OtherDenoisingMain ${ITK_LIBRARIES})

