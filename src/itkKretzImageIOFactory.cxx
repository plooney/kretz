#include "itkKretzImageIOFactory.h"
#include "itkKretzImageIO.h"
#include "itkVersion.h"

namespace itk
{
KretzImageIOFactory::KretzImageIOFactory()
{
  this->RegisterOverride( "itkImageIOBase",
                          "itkKretzImageIO",
                          "Kretz Image IO",
                          1,
                          CreateObjectFunction< KretzImageIO >::New() );
}

KretzImageIOFactory::~KretzImageIOFactory()
{}

const char * KretzImageIOFactory::GetITKSourceVersion() const
{
  return ITK_SOURCE_VERSION;
}

const char * KretzImageIOFactory::GetDescription() const
{
  return "Kretz ImageIO Factory, allows the loading of DICOM images into Insight";
}

// Undocumented API used to register during static initialization.
// DO NOT CALL DIRECTLY.

static bool KretzImageIOFactoryHasBeenRegistered;

void ITKIOKretz_EXPORT KretzImageIOFactoryRegister__Private(void)
{
  if( ! KretzImageIOFactoryHasBeenRegistered )
    {
    KretzImageIOFactoryHasBeenRegistered = true;
    KretzImageIOFactory::RegisterOneFactory();
    }
}

} // end namespace itk
