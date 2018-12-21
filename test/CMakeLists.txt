itk_module_test()

set(IOKretzTests
  itkToroidalTransformTest.cxx
  )

CreateTestDriver(IOKretz "${IOKretz-Test_LIBRARIES}" "${IOKretzTests}")

itk_add_test(NAME itkToroidalTransformTest
  COMMAND IOKretzTestDriver
    #--compare
    #${ITK_TEST_OUTPUT_DIR}/itkMyFilterTestOutput.mha
    #DATA{Baseline/itkMyFilterTestOutput.mha}
  itkToroidalTransformTest
    ${ITK_TEST_OUTPUT_DIR}/testImageCartesian.nii.gz
    ${ITK_TEST_OUTPUT_DIR}/testImageToroidal.nii.gz
    ${ITK_TEST_OUTPUT_DIR}/testImageCartesian2.nii.gz
  )
