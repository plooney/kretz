set(3DUltrasound_SRCS
	itkKretzImageIO.cxx 
	itkKretzImageIOFactory.cxx
	)

itk_module_add_library(ITK3DUltrasound ${3DUltrasound_SRCS})
