cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)
cmake_policy(VERSION 3.5.1)

project(ITK3DUltrasound)

set(ToroidalTransform_LIBRARIES ToroidalTransform)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )
set(CMAKE_CXX_STANDARD 14)
if(USE_CLANG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -Wno-inconsistent-missing-override -Wno-missing-braces")
endif(USE_CLANG)

if(NOT ITK_SOURCE_DIR)
	find_package(ITK REQUIRED)
	list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
	#include(ITKModuleExternal)

	enable_testing()
	include(Dart)

	find_package( Boost 1.58 COMPONENTS program_options REQUIRED )
	find_package(ITK REQUIRED)
	include_directories(src include ${Boost_INCLUDE_DIR} )

	set(Libraries ${ITK_LIBRARIES})

	set(SOURCE
		src/itkKretzWriter.cxx
		src/itkKretzImageIO.cxx
		#src/itkCartesianToToroidalTransform.cxx
		#src/itkToroidalToCartesianTransform.cxx
		)
	include( ${ITK_USE_FILE} )

	add_executable(KretzFileWriter ${SOURCE} )
	target_link_libraries(KretzFileWriter ${ITK_LIBRARIES} ${Boost_LIBRARIES})

	add_executable(KretzConverter
		src/itkKretzConverter.cxx
		src/itkKretzImageIO.cxx
		#src/itkCartesianToToroidalTransform.cxx
		#src/itkToroidalToCartesianTransform.cxx
		)

	target_link_libraries(KretzConverter ${Libraries} ${Boost_LIBRARIES})


	if(BUILD_TESTING)
		add_executable(ToroidalTransformTest
			#src/itkToroidalToCartesianTransform.cxx
			#src/itkCartesianToToroidalTransform.cxx
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

else()
	itk_module_impl()
endif()

include(ExternalProject)

