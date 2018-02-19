cmake_minimum_required(VERSION 2.8)
INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

PROJECT(ToroidalTransform)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")

ENABLE_TESTING()
INCLUDE(Dart)

FIND_PACKAGE( Boost 1.58 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES(src include ${Boost_INCLUDE_DIR} )

SET(Required_Packages
    ITK
    )
FOREACH(Package ${Required_Packages})
LOADPACKAGE(${Package})
ENDFOREACH(Package)

set(Libraries ${ITK_LIBRARIES})
#add_subdirectory(src)
#add_subdirectory(include)

SET(SOURCE
    src/KretzFileWriter.cxx
    src/itkKretzImageIO.cxx
    )

SET(HEADERS
   include/itkKretzImageIO.h
   )

add_executable(KretzFileWriter ${SOURCE} )
target_link_libraries(KretzFileWriter ${} ${ITK_LIBRARIES} ${Boost_LIBRARIES})

add_executable(KretzConverter
    src/itkKretzConverter.cxx
    src/itkKretzImageIO.cxx
    src/itkCartesianToToroidalTransform.txx
    src/itkToroidalToCartesianTransform.txx
    )

target_link_libraries(KretzConverter ${Libraries} ${Boost_LIBRARIES})

add_executable(ToroidalTransformTest
    src/itkToroidalToCartesianTransform.txx
    src/itkCartesianToToroidalTransform.txx
    src/itkToroidalTransformTest.cxx
    )

target_link_libraries(ToroidalTransformTest ${Libraries})

ADD_TEST(ToroidalTransformTest ToroidalTransformTest ToroidalTransformTest)
