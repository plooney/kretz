/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/

#include "itkVersion.h"
#include "itkKretzImageIO.h"
#include "itkIOCommon.h"
#include "itkArray.h"
#include "itkByteSwapper.h"
#include "vnl/vnl_cross.h"

#include "itkMetaDataObject.h"

#include "itksys/SystemTools.hxx"
#include "itksys/Base64.h"

#include <fstream>
#include <sstream>

namespace itk
{

/*class InternalHeader
{
public:
  InternalHeader():m_Header(ITK_NULLPTR) {}
  ~InternalHeader()
  {
    delete m_Header;
  }
//gdcm::File *m_Header;
};
*/

KretzImageIO::KretzImageIO()
{
  /*this->m_KRETZHeader = new InternalHeader;
  this->SetNumberOfDimensions(3); //needed for getting the 3 coordinates of
                                  // the origin, even if it is a 2D slice.
  m_ByteOrder = LittleEndian;     //default
  m_FileType = Binary;            //default...always true
  m_RescaleSlope = 1.0;
  m_RescaleIntercept = 0.0;
  // UIDPrefix is the ITK root id tacked with a ".1"
  // allowing to designate a subspace of the id space for ITK generated DICOM
  m_UIDPrefix = "1.2.826.0.1.3680043.2.1125." "1";

  // Purely internal use, no user access:
  m_StudyInstanceUID = "";
  m_SeriesInstanceUID = "";
  m_FrameOfReferenceInstanceUID = "";

  m_KeepOriginalUID = false;

  m_LoadPrivateTags = false;


  // by default assume that images will be 2D.
  // This number is updated according the information
  // received through the MetaDataDictionary
  m_GlobalNumberOfDimensions = 2;
  // By default use JPEG2000. For legacy system, one should prefer JPEG since
  // JPEG2000 was only recently added to the DICOM standard
  m_CompressionType = JPEG2000;
  */
  this->SetNumberOfDimensions(3); //needed for getting the 3 coordinates of
  this->SetNumberOfComponents(1);
  m_ByteOrder = LittleEndian;     //default
  m_ComponentType = UCHAR;
  m_PixelType = SCALAR;
  m_GlobalNumberOfDimensions = 3;
}

KretzImageIO::~KretzImageIO()
{
  //delete this->m_KretzHeader;
}

bool KretzImageIO::CanReadFile(const char *filename)
{
  std::ifstream file;
  try
    {
    this->OpenFileForReading( file, filename );
    }
  catch( ExceptionObject & )
    {
    return false;
    }

  //
  // sniff for the KRETZ signature 
  //
  bool kretzsig(false);
  {
    file.seekg(0,std::ios_base::beg);
    if(file.fail() || file.eof())
    {
	    return false;
    }
    char buf[14];
    file.read(buf,14);
    if(file.fail())
    {
	    return false;
    }
    //buf[13] = '\0';
    std::string sig(buf);
    if(sig == "KRETZFILE 1.0")
    {
	    kretzsig = true;
	    return true;
    }
  }
  return false;
}

void KretzImageIO::Read(void *buffer)
{

  // ensure file can be opened for reading, before doing any more work
  std::ifstream inputFileStream;
  this->OpenFileForReading( inputFileStream, m_FileName );

  // Find last position for tag search
  inputFileStream.seekg(5,std::ios::end);
  int endOfTagSearch = inputFileStream.tellg();

  inputFileStream.seekg(16,std::ios::beg); // skips header

  // let any exceptions propagate
  //buff = reinterpret_cast<char *>(&tag), sizeof(tag))
  unsigned short tag_shorts[2];
  unsigned int taglength;
  while((inputFileStream.tellg()<endOfTagSearch-8)&&(inputFileStream.tellg()>-1)){
    inputFileStream.read(reinterpret_cast<char *>(&tag_shorts), sizeof(tag_shorts)); 
    inputFileStream.read(reinterpret_cast<char *>(&taglength), sizeof(taglength)); 
    int searchIterator = inputFileStream.tellg(); //added to allow for tag length/data type mismatch;
    bool val_unread = true;
   
    Tag tag(tag_shorts[0],tag_shorts[1]);
    
    if(tag==ImageTag){
	    
        inputFileStream.read(reinterpret_cast<char *>(buffer), taglength); 
	std::cout << "Found image data len: " << taglength  << std::endl;
    } 
    else {
      inputFileStream.seekg(taglength,std::ios::cur);
    }


    //After we've read the data, need to check we're still in line with the file.
    if((int)inputFileStream.tellg()!=searchIterator+taglength){
      //skip any extra padding required to allow for tag length/data size mismatch
      std::cout << "In unwanted loop" << std::endl;
      inputFileStream.seekg(searchIterator,std::ios::beg);
      inputFileStream.seekg(taglength,std::ios::cur);
    }
  }

  inputFileStream.close();
  /*

  itkAssertInDebugAndIgnoreInReleaseMacro( gdcm::ImageHelper::GetForceRescaleInterceptSlope() );
  gdcm::ImageReader reader;
  reader.SetFileName( m_FileName.c_str() );
  if ( !reader.Read() )
    {
    itkExceptionMacro(<< "Cannot read requested file");
    }

  gdcm::Image & image = reader.GetImage();
#ifndef NDEBUG
  gdcm::PixelFormat pixeltype_debug = image.GetPixelFormat();
  itkAssertInDebugAndIgnoreInReleaseMacro(image.GetNumberOfDimensions() == 2 || image.GetNumberOfDimensions() == 3);
#endif
  SizeValueType len = image.GetBufferLength();

  // I think ITK only allow RGB image by pixel (and not by plane)
  if ( image.GetPlanarConfiguration() == 1 )
    {
    gdcm::ImageChangePlanarConfiguration icpc;
    icpc.SetInput(image);
    icpc.SetPlanarConfiguration(0);
    icpc.Change();
    image = icpc.GetOutput();
    }

  gdcm::PhotometricInterpretation pi = image.GetPhotometricInterpretation();
  if ( pi == gdcm::PhotometricInterpretation::PALETTE_COLOR )
    {
    gdcm::ImageApplyLookupTable ialut;
    ialut.SetInput(image);
    ialut.Apply();
    image = ialut.GetOutput();
    len *= 3;
    }

  if ( !image.GetBuffer( (char*)pointer ) )
    {
    itkExceptionMacro(<< "Failed to get the buffer!");
    }

  const gdcm::PixelFormat & pixeltype = image.GetPixelFormat();
#ifndef NDEBUG
  // ImageApplyLookupTable is meant to change the pixel type for PALETTE_COLOR images
  // (from single values to triple values per pixel)
  if ( pi != gdcm::PhotometricInterpretation::PALETTE_COLOR )
    {
    itkAssertInDebugAndIgnoreInReleaseMacro( pixeltype_debug == pixeltype );
    }
#endif

  if ( m_RescaleSlope != 1.0 || m_RescaleIntercept != 0.0 )
    {
    gdcm::Rescaler r;
    r.SetIntercept(m_RescaleIntercept);
    r.SetSlope(m_RescaleSlope);
    r.SetPixelFormat(pixeltype);
    gdcm::PixelFormat outputpt = r.ComputeInterceptSlopePixelType();
    char *            copy = new char[len];
    memcpy(copy, (char *)pointer, len);
    r.Rescale( (char *)pointer, copy, len );
    delete[] copy;
    // WARNING: sizeof(Real World Value) != sizeof(Stored Pixel)
    len = len * outputpt.GetPixelSize() / pixeltype.GetPixelSize();
    }

#ifndef NDEBUG
  // \postcondition
  // Now that len was updated (after unpacker 12bits -> 16bits, rescale...) ,
  // can now check compat:
  const SizeValueType numberOfBytesToBeRead =
    static_cast< SizeValueType >( this->GetImageSizeInBytes() );
  itkAssertInDebugAndIgnoreInReleaseMacro(numberOfBytesToBeRead == len);   // programmer error
#endif
  */
}

void KretzImageIO::ReadImageInformation()
{
  std::ifstream inputFileStream;
  this->OpenFileForReading( inputFileStream, m_FileName );

  // Find last position for tag search
  inputFileStream.seekg(5,std::ios::end);
  int endOfTagSearch = inputFileStream.tellg();

  inputFileStream.seekg(16,std::ios::beg); // skips header

  // let any exceptions propagate
  //buff = reinterpret_cast<char *>(&tag), sizeof(tag))
  unsigned short tag_shorts[2];
  unsigned int taglength;
  while((inputFileStream.tellg()<endOfTagSearch-8)&&(inputFileStream.tellg()>-1)){
    inputFileStream.read(reinterpret_cast<char *>(&tag_shorts), sizeof(tag_shorts)); 
    inputFileStream.read(reinterpret_cast<char *>(&taglength), sizeof(taglength)); 
    int searchIterator = inputFileStream.tellg(); //added to allow for tag length/data type mismatch;
    bool val_unread = true;
   
    Tag tag(tag_shorts[0],tag_shorts[1]);
    
    if(tag==PatientTag){
	char patient_name_c_array[taglength]; 
        inputFileStream.read(reinterpret_cast<char *>(patient_name_c_array), taglength); 
	std::cout << std::string(patient_name_c_array) << std::endl;
    }
    else if(tag==DimensionXTag){
        unsigned short dimension;
        inputFileStream.read(reinterpret_cast<char *>(&dimension), taglength); 
	std::cout << dimension << std::endl;
	this->SetDimensions(0,dimension);
    } 
    else if(tag==DimensionYTag){
        unsigned short dimension;
        inputFileStream.read(reinterpret_cast<char *>(&dimension), taglength); 
	std::cout << dimension << std::endl;
	this->SetDimensions(1,dimension);
    } 
    else if(tag==DimensionZTag){
        unsigned short dimension;
        inputFileStream.read(reinterpret_cast<char *>(&dimension), taglength); 
	std::cout << dimension << std::endl;
	this->SetDimensions(2,dimension);
    } 
    else if(tag==ResolutionTag){
        double value;
        inputFileStream.read(reinterpret_cast<char *>(&value), taglength); 
	this->m_Resolution = value *1000.0;
    } 
    else if(tag==Offset1Tag){
        double value;
        inputFileStream.read(reinterpret_cast<char *>(&value), taglength); 
	this->m_Offset1 = value;
    } 
    else if(tag==Offset2Tag){
        double value;
        inputFileStream.read(reinterpret_cast<char *>(&value), taglength); 
	this->m_Offset2 = value;
    } 
    else if(tag==Angles1Tag){
        int len = taglength/sizeof(double);
	double * angles = new double[len];
        inputFileStream.read(reinterpret_cast<char *>(angles), sizeof( double ) * len); 
        double amin = angles[0];
        double amax = angles[len-1];
	double bCentre = (amin+amax)/2;
	for(int i=0; i<len;i++){
            double angle = angles[i]-bCentre;
	    this->m_TableAngles2.push_back(std::make_pair(angle, i));
	}
	delete angles;
    } 
    else if(tag==Angles2Tag){
        int len = taglength/sizeof(double);
	double * angles = new double[len];
        inputFileStream.read(reinterpret_cast<char *>(angles), sizeof( double ) * len); 
        double amin = angles[0];
        double amax = angles[len-1];
	double bCentre = (amin+amax)/2;
	for(int i=0; i<len;i++){
            double angle = angles[i]-bCentre;
	    this->m_TableAngles1.push_back(std::make_pair(angle, i));
	}
	delete angles;
    } 
    else {
      inputFileStream.seekg(taglength,std::ios::cur);
    }


    //After we've read the data, need to check we're still in line with the file.
    if((int)inputFileStream.tellg()!=searchIterator+taglength){
      //skip any extra padding required to allow for tag length/data size mismatch
      std::cout << "In unwanted loop" << std::endl;
      inputFileStream.seekg(searchIterator,std::ios::beg);
      inputFileStream.seekg(taglength,std::ios::cur);
    }
  }
  m_rBstart = m_Offset1*m_Resolution;
  m_rD = -m_Offset2*m_Resolution;

  inputFileStream.close();
  MetaDataDictionary & dico = this->GetMetaDataDictionary();

}

/*
void KretzImageIO::InternalReadImageInformation()
{
  // ensure file can be opened for reading, before doing any more work
  // In general this should be relatively safe to assume
  gdcm::ImageHelper::SetForceRescaleInterceptSlope(true);

  gdcm::ImageReader reader;
  reader.SetFileName( m_FileName.c_str() );
  if ( !reader.Read() )
    {
    itkExceptionMacro(<< "Cannot read requested file");
    }
  const gdcm::Image &   image = reader.GetImage();
  const gdcm::File &    f = reader.GetFile();
  const gdcm::DataSet & ds = f.GetDataSet();
  const unsigned int *  dims = image.GetDimensions();

  const gdcm::PixelFormat & pixeltype = image.GetPixelFormat();
  switch ( pixeltype )
    {
    case gdcm::PixelFormat::INT8:
      m_InternalComponentType = ImageIOBase::CHAR; // Is it signed char ?
      break;
    case gdcm::PixelFormat::UINT8:
      m_InternalComponentType = ImageIOBase::UCHAR;
      break;
    case gdcm::PixelFormat::INT12:
      m_InternalComponentType = ImageIOBase::SHORT;
      break;
    case gdcm::PixelFormat::UINT12:
      m_InternalComponentType = ImageIOBase::USHORT;
      break;
    case gdcm::PixelFormat::INT16:
      m_InternalComponentType = ImageIOBase::SHORT;
      break;
    case gdcm::PixelFormat::UINT16:
      m_InternalComponentType = ImageIOBase::USHORT;
      break;
    // RT / SC have 32bits
    case gdcm::PixelFormat::INT32:
      m_InternalComponentType = ImageIOBase::INT;
      break;
    case gdcm::PixelFormat::UINT32:
      m_InternalComponentType = ImageIOBase::UINT;
      break;
    //case gdcm::PixelFormat::FLOAT16: // TODO
    case gdcm::PixelFormat::FLOAT32:
      m_InternalComponentType = ImageIOBase::FLOAT;
      break;
    case gdcm::PixelFormat::FLOAT64:
      m_InternalComponentType = ImageIOBase::DOUBLE;
      break;
    default:
      itkExceptionMacro("Unhandled PixelFormat: " << pixeltype);
    }

  m_RescaleIntercept = image.GetIntercept();
  m_RescaleSlope = image.GetSlope();
  gdcm::Rescaler r;
  r.SetIntercept(m_RescaleIntercept);
  r.SetSlope(m_RescaleSlope);
  r.SetPixelFormat(pixeltype);
  gdcm::PixelFormat::ScalarType outputpt = r.ComputeInterceptSlopePixelType();

  itkAssertInDebugAndIgnoreInReleaseMacro(pixeltype <= outputpt);

  m_ComponentType = UNKNOWNCOMPONENTTYPE;
  switch ( outputpt )
    {
    case gdcm::PixelFormat::INT8:
      m_ComponentType = ImageIOBase::CHAR; // Is it signed char ?
      break;
    case gdcm::PixelFormat::UINT8:
      m_ComponentType = ImageIOBase::UCHAR;
      break;
    case gdcm::PixelFormat::INT12:
      m_ComponentType = ImageIOBase::SHORT;
      break;
    case gdcm::PixelFormat::UINT12:
      m_ComponentType = ImageIOBase::USHORT;
      break;
    case gdcm::PixelFormat::INT16:
      m_ComponentType = ImageIOBase::SHORT;
      break;
    case gdcm::PixelFormat::UINT16:
      m_ComponentType = ImageIOBase::USHORT;
      break;
    // RT / SC have 32bits
    case gdcm::PixelFormat::INT32:
      m_ComponentType = ImageIOBase::INT;
      break;
    case gdcm::PixelFormat::UINT32:
      m_ComponentType = ImageIOBase::UINT;
      break;
    //case gdcm::PixelFormat::FLOAT16: // TODO
    case gdcm::PixelFormat::FLOAT32:
      m_ComponentType = ImageIOBase::FLOAT;
      break;
    case gdcm::PixelFormat::FLOAT64:
      m_ComponentType = ImageIOBase::DOUBLE;
      break;
    default:
      itkExceptionMacro("Unhandled PixelFormat: " << outputpt);
    }

  m_NumberOfComponents = pixeltype.GetSamplesPerPixel();
  if ( image.GetPhotometricInterpretation() ==
       gdcm::PhotometricInterpretation::PALETTE_COLOR )
    {
    itkAssertInDebugAndIgnoreInReleaseMacro(m_NumberOfComponents == 1);
    // TODO: need to do the LUT ourself...
    //itkExceptionMacro(<< "PALETTE_COLOR is not implemented yet");
    // AFAIK ITK user don't care about the palette so always apply it and fake a
    // RGB image for them
    m_NumberOfComponents = 3;
    }
  if ( m_NumberOfComponents == 1 )
    {
    this->SetPixelType(SCALAR);
    }
  else
    {
    this->SetPixelType(RGB); // What if image is YBR ? This is a problem since
                             // integer conversion is lossy
    }

  // set values in case we don't find them
  //this->SetNumberOfDimensions(  image.GetNumberOfDimensions() );
  m_Dimensions[0] = dims[0];
  m_Dimensions[1] = dims[1];

  double spacing[3];

  //
  //
  // This is a WORKAROUND for a bug in GDCM -- in
  // ImageHeplper::GetSpacingTagFromMediaStorage it was not
  // handling some MediaStorage types
  // so we have to punt here.
  gdcm::MediaStorage ms;

  ms.SetFromFile(f);
  switch(ms)
    {
    case gdcm::MediaStorage::HardcopyGrayscaleImageStorage:
    case gdcm::MediaStorage::GEPrivate3DModelStorage:
    case gdcm::MediaStorage::Philips3D:
    case gdcm::MediaStorage::VideoEndoscopicImageStorage:
    case gdcm::MediaStorage::UltrasoundMultiFrameImageStorage:
    case gdcm::MediaStorage::UltrasoundImageStorage: // ??
    case gdcm::MediaStorage::UltrasoundImageStorageRetired:
    case gdcm::MediaStorage::UltrasoundMultiFrameImageStorageRetired:
      {
      std::vector<double> sp;
      gdcm::Tag spacingTag(0x0028,0x0030);
      if(ds.FindDataElement(spacingTag) &&
         !ds.GetDataElement(spacingTag).IsEmpty())
        {
        gdcm::DataElement de = ds.GetDataElement(spacingTag);
        const gdcm::Global &g = gdcm::GlobalInstance;
        const gdcm::Dicts &dicts = g.GetDicts();
        const gdcm::DictEntry &entry = dicts.GetDictEntry(de.GetTag());
        const gdcm::VR & vr = entry.GetVR();
        switch(vr)
          {
          case gdcm::VR::DS:
            {
            std::stringstream m_Ss;

            gdcm::Element<gdcm::VR::DS,gdcm::VM::VM1_n> m_El;

            const gdcm::ByteValue *                     bv = de.GetByteValue();
            assert( bv );

            std::string s = std::string( bv->GetPointer(), bv->GetLength() );

            m_Ss.str( s );
            // Stupid file: CT-MONO2-8-abdo.dcm
            // The spacing is something like that: [0.2\0\0.200000]
            // I would need to throw an expection that VM is not compatible
            m_El.SetLength( entry.GetVM().GetLength() * entry.GetVR().GetSizeof() );
            m_El.Read( m_Ss );

            assert( m_El.GetLength() == 2 );
            for(unsigned long i = 0; i < m_El.GetLength(); ++i)
              sp.push_back( m_El.GetValue(i) );
            std::swap( sp[0], sp[1]);
            assert( sp.size() == (unsigned int)entry.GetVM() );
            }
            break;
          case gdcm::VR::IS:
            {
            std::stringstream m_Ss;

            gdcm::Element<gdcm::VR::IS,gdcm::VM::VM1_n> m_El;

            const gdcm::ByteValue *bv = de.GetByteValue();
            assert( bv );

            std::string s = std::string( bv->GetPointer(), bv->GetLength() );
            m_Ss.str( s );
            m_El.SetLength( entry.GetVM().GetLength() * entry.GetVR().GetSizeof() );
            m_El.Read( m_Ss );
            for(unsigned long i = 0; i < m_El.GetLength(); ++i)
              sp.push_back( m_El.GetValue(i) );
            assert( sp.size() == (unsigned int)entry.GetVM() );
            }
            break;
          default:
            assert(0);
            break;
          }
        spacing[0] = sp[0];
        spacing[1] = sp[1];
        }
      else
        {
        spacing[0] = 1.0;
        spacing[1] = 1.0;
        }
      spacing[2] = 1.0; // punt?
      }
      break;
    default:
      {
      const double *sp;
      sp = image.GetSpacing();
      spacing[0] = sp[0];
      spacing[1] = sp[1];
      spacing[2] = sp[2];
      }
    break;
    }

  const double *origin = image.GetOrigin();
  for(unsigned i = 0; i < 3; ++i)
    {
    m_Spacing[i] = spacing[i];
    m_Origin[i] = origin[i];
    }
  if ( image.GetNumberOfDimensions() == 3 )
    {
    m_Dimensions[2] = dims[2];
    }
  else
    {
    m_Dimensions[2] = 1;
    }

  const double *       dircos = image.GetDirectionCosines();
  vnl_vector< double > rowDirection(3), columnDirection(3);
  rowDirection[0] = dircos[0];
  rowDirection[1] = dircos[1];
  rowDirection[2] = dircos[2];

  columnDirection[0] = dircos[3];
  columnDirection[1] = dircos[4];
  columnDirection[2] = dircos[5];

  vnl_vector< double > sliceDirection = vnl_cross_3d(rowDirection, columnDirection);

  // orthogonalize
  sliceDirection.normalize();
  rowDirection = vnl_cross_3d(columnDirection,sliceDirection).normalize();
  columnDirection = vnl_cross_3d(sliceDirection,rowDirection);

  this->SetDirection(0, rowDirection);
  this->SetDirection(1, columnDirection);
  this->SetDirection(2, sliceDirection);

  //Now copying the gdcm dictionary to the itk dictionary:
  MetaDataDictionary & dico = this->GetMetaDataDictionary();

  // in the case of re-use by ImageSeriesReader, clear the dictionary
  // before populating it.
  dico.Clear();

  gdcm::StringFilter sf;
  sf.SetFile(f);
  gdcm::DataSet::ConstIterator it = ds.Begin();

  // Copy of the header->content
  for (; it != ds.End(); ++it )
    {
    const gdcm::DataElement & ref = *it;
    const gdcm::Tag &         tag = ref.GetTag();
    // Compute VR from the toplevel file, and the currently processed dataset:
    gdcm::VR vr = gdcm::DataSetHelper::ComputeVR(f, ds, tag);

    // Process binary field and encode them as mime64: only when we do not know
    // of any better
    // representation. VR::US is binary, but user want ASCII representation.
    if ( vr & ( gdcm::VR::OB | gdcm::VR::OF | gdcm::VR::OW | gdcm::VR::SQ | gdcm::VR::UN ) )
      {
      // itkAssertInDebugAndIgnoreInReleaseMacro( vr & gdcm::VR::VRBINARY );
      if ( (m_LoadPrivateTags || tag.IsPublic()) && vr != gdcm::VR::SQ
           && tag != gdcm::Tag(0x7fe0, 0x0010) // && vr != gdcm::VR::UN )
        {
        const gdcm::ByteValue *bv = ref.GetByteValue();
        if ( bv )
          {
          // base64 streams have to be a multiple of 4 bytes in length
          int encodedLengthEstimate = 2 * bv->GetLength();
          encodedLengthEstimate = ( ( encodedLengthEstimate / 4 ) + 1 ) * 4;

          char *       bin = new char[encodedLengthEstimate];
          unsigned int encodedLengthActual = static_cast< unsigned int >(
            itksysBase64_Encode(
              (const unsigned char *)bv->GetPointer(),
              static_cast< SizeValueType >( bv->GetLength() ),
              (unsigned char *)bin,
              static_cast< int >( 0 ) ) );
          std::string encodedValue(bin, encodedLengthActual);
          EncapsulateMetaData< std::string >(dico, tag.PrintAsPipeSeparatedString(), encodedValue);
          delete[] bin;
          }
        }
      }
    else // if ( vr & gdcm::VR::VRASCII ) 
      {
      // Only copying field from the public DICOM dictionary
      if ( m_LoadPrivateTags || tag.IsPublic() )
        {
        EncapsulateMetaData< std::string >( dico, tag.PrintAsPipeSeparatedString(), sf.ToString(tag) );
        }
      }
    }
}
*/

bool KretzImageIO::CanWriteFile(const char *name)
{
  return false;
}

void KretzImageIO::WriteImageInformation()
{}

void KretzImageIO::Write(const void *buffer)
{
    {
    itkExceptionMacro(<< "Kretz file writing is not supported.");
    }
}
/*
// Convenience methods to query patient and scanner information. These
// methods are here for compatibility with the DICOMImageIO2 class.
void KretzImageIO::GetPatientName(char *name)
{
  MetaDataDictionary & dict = this->GetMetaDataDictionary();

  ExposeMetaData< std::string >(dict, "0010|0010", m_PatientName);
  strcpy ( name, m_PatientName.c_str() );
}


bool KretzImageIO::GetValueFromTag(const std::string & tag, std::string & value)
{
  MetaDataDictionary & dict = this->GetMetaDataDictionary();

  std::string tag_lower = tag;
  std::transform( tag_lower.begin(), tag_lower.end(), tag_lower.begin(),
                  static_cast<int(*)(int)>( ::tolower ) );

  return ExposeMetaData< std::string >(dict, tag_lower, value);
}

bool KretzImageIO::GetLabelFromTag(const std::string & tag,
                                  std::string & labelId)
{
  gdcm::Tag t;
  if ( t.ReadFromPipeSeparatedString( tag.c_str() ) && t.IsPublic() )
    {
    const gdcm::Global &    g = gdcm::Global::GetInstance();
    const gdcm::Dicts &     dicts = g.GetDicts();
    const gdcm::DictEntry & entry = dicts.GetDictEntry(t);
    labelId = entry.GetName();
    return true;
    }
  return false;
}
*/

void KretzImageIO::PrintSelf(std::ostream & os, Indent indent) const
{
/*
  Superclass::PrintSelf(os, indent);
  os << indent << "Internal Component Type: " << this->GetComponentTypeAsString(m_InternalComponentType)
     << std::endl;
  os << indent << "RescaleSlope: " << m_RescaleSlope << std::endl;
  os << indent << "RescaleIntercept: " << m_RescaleIntercept << std::endl;
  os << indent << "KeepOriginalUID:" << ( m_KeepOriginalUID ? "On" : "Off" ) << std::endl;
  os << indent << "LoadPrivateTags:" << ( m_LoadPrivateTags ? "On" : "Off" ) << std::endl;
  os << indent << "UIDPrefix: " << m_UIDPrefix << std::endl;
  os << indent << "StudyInstanceUID: " << m_StudyInstanceUID << std::endl;
  os << indent << "SeriesInstanceUID: " << m_SeriesInstanceUID << std::endl;
  os << indent << "FrameOfReferenceInstanceUID: " << m_FrameOfReferenceInstanceUID << std::endl;
  os << indent << "CompressionType:" << m_CompressionType << std::endl;
*/
}
} // end namespace itk
