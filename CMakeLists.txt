cmake_minimum_required(VERSION 2.8)
include(${CMAKE_SOURCE_DIR}/IJMacros.txt)
include(ExternalProject)

project(ToroidalTransform)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )
set(CMAKE_CXX_STANDARD 14)

enable_testing()
include(Dart)

find_package( Boost 1.58 COMPONENTS program_options REQUIRED )
include_directories(src include ${Boost_INCLUDE_DIR} )

SET(Required_Packages
	ITK
	)
foreach(Package ${Required_Packages})
	loadpackage(${Package})
endforeach(Package)

set(Libraries ${ITK_LIBRARIES})
#add_subdirectory(src)
#add_subdirectory(include)

set(SOURCE
	src/itkKretzWriter.cxx
	src/itkKretzImageIO.cxx
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

if(BUILD_TESTING)
	add_executable(ToroidalTransformTest
		src/itkToroidalToCartesianTransform.txx
		src/itkCartesianToToroidalTransform.txx
		test/itkToroidalTransformTest.cxx
		)

	target_link_libraries(ToroidalTransformTest ${Libraries})

	# Test writing the toroidal volume out to a nii.gz file format, geometry of the file is determined by angles in  
	add_test(KretzFileTest KretzFileWriter -i test/test.vol -o test/test_toroidal.nii.gz) 

	# Test conversion of toroidal Kretz fie to cartesian coordinates
	add_test(KretzConverterTest KretzConverter -i test/test.vol -o test/test_cartesian.nii.gz -r 0.6 0.6 0.6) 

	# Test conversion of cartesian file to toroidal nii.gz with geometry specified in the Kretz file
	add_test(KretzFileTestCartesian KretzFileWriter -i test/test.vol -o test/test_toroidal2.nii.gz -c test/test_cartesian.nii.gz) 

	# Test conversion of created toroidal volume to cartesian and back again
	add_test(ToroidalTransformTest ToroidalTransformTest)

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

find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API        
documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
	if(NOT DOXYGEN_FOUND)
		message(FATAL_ERROR "Doxygen is needed to build the documentation.")
	endif()

	set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/kretz.doxygen)
	set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)

	configure_file(${doxyfile_in} ${doxyfile} @ONLY)

	add_custom_target(doc
		${DOXYGEN_EXECUTABLE} doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
		)
endif()
