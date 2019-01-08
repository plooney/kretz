cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)
cmake_policy(VERSION 3.10.2)
set(CMAKE_CXX_STANDARD 14)

project(ITK3DUltrasound)

find_package(ITK REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})

include(CTest)

find_package(Boost 1.58 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
find_package(ITK REQUIRED
  COMPONENTS
    IOKretz
    ITKImageIO
    ITKImageGrid
    ITKImageIntensity

  )
include(${ITK_USE_FILE})

add_executable(KretzFileWriter itkKretzWriter.cxx)
target_link_libraries(KretzFileWriter ${ITK_LIBRARIES} ${Boost_LIBRARIES})

add_executable(KretzConverter itkKretzConverter.cxx)
target_link_libraries(KretzConverter ${ITK_LIBRARIES} ${Boost_LIBRARIES})


if(BUILD_TESTING)
  # Test writing the toroidal volume out to a nii.gz file format, geometry of the file is determined by angles in
  add_test(KretzFileTest KretzFileWriter -i test/test.vol -o test/test_toroidal.nii.gz)

  # Test conversion of toroidal Kretz fie to cartesian coordinates
  add_test(KretzConverterTest KretzConverter -i test/test.vol -o test/test_cartesian.nii.gz -r 0.6 0.6 0.6)

  # Test conversion of cartesian file to toroidal nii.gz with geometry specified in the Kretz file
  add_test(KretzFileTestCartesian KretzFileWriter -i test/test.vol -o test/test_toroidal2.nii.gz -c test/test_cartesian.nii.gz)

  file(DOWNLOAD https://data.kitware.com/api/v1/item/5b5712cd8d777f06857bfd77/download test/test.vol
          EXPECTED_MD5 32631727d0202dc825f12418806f2f13
          )
  file(DOWNLOAD https://data.kitware.com/api/v1/item/5b5712e18d777f06857bfd7a/download expected_test_output/test_cartesian.nii.gz
          EXPECTED_MD5 a26bdeec4e43a122b68c7bbfcd35d54c
          )
  file(DOWNLOAD https://data.kitware.com/api/v1/item/5b5712e88d777f06857bfd7d/download expected_test_output/test_toroidal.nii.gz
          EXPECTED_MD5 1a307bbcd05ad9d5bbd02ef2fbdeaf50
          )
endif()
