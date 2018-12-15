set(IOKretz_SRCS
  itkKretzImageIO.cxx
  itkKretzImageIOFactory.cxx
  )

itk_module_add_library(IOKretz ${IOKretz_SRCS})
