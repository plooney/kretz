# the top-level README is used for describing this module, just
# re-used it for documentation here
get_filename_component(MY_CURRENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(READ "${MY_CURRENT_DIR}/README.md" DOCUMENTATION)

itk_module(IOKretz
  DEPENDS
    ITKCommon
    ITKIOImageBase
  COMPILE_DEPENDS
    ITKTransform
    ITKImageFunction
    ITKMesh
  TEST_DEPENDS
    ITKTestKernel
    ITKImageGrid
    ITKImageStatistics
  FACTORY_NAMES
    ImageIO::Kretz
  DESCRIPTION
    "${DOCUMENTATION}"
  EXCLUDE_FROM_DEFAULT
  ENABLE_SHARED
)
