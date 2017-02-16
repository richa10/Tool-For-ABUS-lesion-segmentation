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


# add_executable(SegmentationABUS_fixedParam SegmentationABUS_fixedParam.cxx )
# target_link_libraries(SegmentationABUS_fixedParam ${ITK_LIBRARIES})

# add_executable(Segmentation_ABUS_Itera_sigma Segmentation_ABUS_Itera_sigma.cxx )
# target_link_libraries(Segmentation_ABUS_Itera_sigma ${ITK_LIBRARIES})

# add_executable(imagediff imagediff.cxx )
# target_link_libraries(imagediff ${ITK_LIBRARIES})


# add_executable(gradientcnr gradientcnr.cxx )
# target_link_libraries(gradientcnr ${ITK_LIBRARIES})

add_executable(AnisotropicMain AnisotropicMain.cxx )
target_link_libraries(AnisotropicMain ${ITK_LIBRARIES})

add_executable(ADLBRMain ADLBRMain.cxx )
target_link_libraries(ADLBRMain ${ITK_LIBRARIES})

# add_executable(AnisotropicMainTao AnisotropicMainTao.cxx )
# target_link_libraries(AnisotropicMainTao ${ITK_LIBRARIES})

# add_executable(ADLBRMainTao ADLBRMainTao.cxx )
# target_link_libraries(ADLBRMainTao ${ITK_LIBRARIES})

add_executable(FileSize FileSize.cxx )
target_link_libraries(FileSize ${ITK_LIBRARIES})

add_executable(ADLBRMainOld1411 ADLBRMainOld1411.cxx )
target_link_libraries(ADLBRMainOld1411 ${ITK_LIBRARIES})

add_executable(ADGaussian ADGaussian.cxx )
target_link_libraries(ADGaussian ${ITK_LIBRARIES})
