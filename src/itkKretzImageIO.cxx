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

// anonymous namespace
namespace
{

/**
 *Class to abstract a group and element for reading a Kretz file.
 */
class Tag
{
public:
  Tag();
  Tag(unsigned short group, unsigned short element){
    this->group = group;
    this->element = element;
  };
  ~Tag(){};
  unsigned short group;
  unsigned short element;
  bool operator==(Tag& other) const
  {
    if(this == &other) return true; //This is the pointer for
    if(this->group == other.group && this->element==other.element) return true;
    else return false;
  }
  operator std::string() const {
    std::stringstream stream;
    stream << std::hex << this->group;
    std::string result1( stream.str() );

    stream << std::hex << this->element;
    std::string result2( stream.str() );
    std::string result(result1 + "|" + result2);

    return result;
  }
};

inline
std::ostream & operator<<(std::ostream & Str, Tag const & v) {
    std::string string(v);
    Str << string;
    return Str;
};
static Tag PatientTag(0x0110, 0x0002);
static Tag CineFramesTag(0xd400, 0x0001);
static Tag SizeFramesTag(0xd400, 0x0002);
static Tag TimingFramesTag(0xd400, 0x0005);

static Tag DimensionXTag(0xc000, 0x0001);
static Tag DimensionYTag(0xc000, 0x0002);
static Tag DimensionZTag(0xc000, 0x0003);
static Tag ResolutionTag(0xc100, 0x0001);
static Tag Offset1Tag(0xc200, 0x0001);
static Tag Offset2Tag(0xc200, 0x0002);
static Tag AnglesPhiTag(0xc300, 0x0001);
static Tag AnglesThetaTag(0xc300, 0x0002);
static Tag ImageTag(0xd000, 0x0001);
static Tag Image4dTag(0xd600, 0x0001);

/*
static Tag PatientTag(0x0110, 0x0002);
static Tag DimensionXTag(0xc000, 0x0201);
static Tag DimensionYTag(0xc000, 0x0202);
static Tag DimensionZTag(0xc000, 0x0203);
static Tag ResolutionTag(0xc100, 0x0201);
static Tag Offset1Tag(0xc200, 0x0201);
static Tag Offset2Tag(0xc200, 0x0202);
static Tag AnglesPhiTag(0xc300, 0x0201);
static Tag AnglesThetaTag(0xc300, 0x0202);
static Tag ImageTag(0xd000, 0x0201);
static Tag CineFramesTag(0xd400, 0x0001);
static Tag SizeFramesTag(0xd400, 0x0002);
static Tag TimingFramesTag(0xd400, 0x0005);
static Tag Image4dTag(0xd600, 0x0201);
*/

static Tag DimensionXTagDoppler(0xc000, 0x0201); static Tag
DimensionYTagDoppler(0xc000, 0x0202); static Tag DimensionZTagDoppler(0xc000,
  0x0203); static Tag ResolutionTagDoppler(0xc100, 0x0201); static Tag
Offset1TagDoppler(0xc200, 0x0201); static Tag Offset2TagDoppler(0xc200,
  0x0202); static Tag AnglesPhiTagDoppler(0xc300, 0x0201); static Tag
AnglesThetaTagDoppler(0xc300, 0x0202); static Tag ImageTagDoppler(0xd000,
  0x0201); static Tag Image4dTagDoppler(0xd600, 0x0201); class InternalHeader;
}

namespace itk
{

KretzImageIO::KretzImageIO()
{
  this->SetNumberOfDimensions(3); //needed for getting the 3 coordinates of
  this->SetNumberOfComponents(1);
  m_ByteOrder = LittleEndian;     //default
  m_ComponentType = UCHAR;
  m_PixelType = SCALAR;
}

KretzImageIO::~KretzImageIO()
{
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

  // sniff for the KRETZ signature 
  {
    file.seekg(0,std::ios_base::beg);
    if(file.fail() || file.eof())
    {
      return false;
    }
    char buf[16];
    file.read(buf,16);
    if(file.fail())
    {
      return false;
    }
    std::string sig(buf);
    if(sig == "KRETZFILE 1.0   ")
    {
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
  unsigned short tag_shorts[2];
  unsigned int taglength;
  while((inputFileStream.tellg()<endOfTagSearch-8)&&(inputFileStream.tellg()>-1))
  {
    inputFileStream.read(reinterpret_cast<char *>(&tag_shorts), sizeof(tag_shorts)); 
    inputFileStream.read(reinterpret_cast<char *>(&taglength), sizeof(taglength)); 

    Tag tag(tag_shorts[0],tag_shorts[1]);

    if(tag==ImageTag && !m_IsDoppler)
    {

      inputFileStream.read(reinterpret_cast<char *>(buffer), taglength); 

    } 
    else if(tag==ImageTagDoppler && m_IsDoppler)
    {

      inputFileStream.read(reinterpret_cast<char *>(buffer), taglength); 

    } 
    else {
      inputFileStream.seekg(taglength,std::ios::cur);
    }

  }

  inputFileStream.close();
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
  MetaDataDictionary & dico = this->GetMetaDataDictionary();

  unsigned short tag_shorts[2];
  unsigned int taglength;
  double offset1, offset2;
  while((inputFileStream.tellg()<endOfTagSearch-8)&&(inputFileStream.tellg()>-1))
  {
    inputFileStream.read(reinterpret_cast<char *>(&tag_shorts), sizeof(tag_shorts)); 
    inputFileStream.read(reinterpret_cast<char *>(&taglength), sizeof(taglength)); 

    Tag tag(tag_shorts[0],tag_shorts[1]);

    if(tag==PatientTag)
    {
      char patient_name_c_array[taglength]; 
      inputFileStream.read(reinterpret_cast<char *>(patient_name_c_array), taglength); 
      std::string tagString(tag);
      EncapsulateMetaData< std::string >( dico, tagString, std::string(patient_name_c_array) );
    }
    else if(tag==DimensionXTag && !m_IsDoppler)
    {
      unsigned short dimension;
      inputFileStream.read(reinterpret_cast<char *>(&dimension), taglength); 
      this->SetDimensions(0,dimension);
    } 
    else if(tag==DimensionYTag && !m_IsDoppler)
    {
      unsigned short dimension;
      inputFileStream.read(reinterpret_cast<char *>(&dimension), taglength); 
      this->SetDimensions(1,dimension);
    } 
    else if(tag==DimensionZTag && !m_IsDoppler)
    {
      unsigned short dimension;
      inputFileStream.read(reinterpret_cast<char *>(&dimension), taglength); 
      this->SetDimensions(2,dimension);
    } 
    else if(tag==ResolutionTag && !m_IsDoppler)
    {
      double value;
      inputFileStream.read(reinterpret_cast<char *>(&value), taglength); 
      this->m_Resolution = value *1000.0;
    } 
    else if(tag==Offset1Tag && !m_IsDoppler)
    {
      double value;
      inputFileStream.read(reinterpret_cast<char *>(&value), taglength); 
      offset1 = value;
    } 
    else if(tag==Offset2Tag && !m_IsDoppler)
    {
      double value;
      inputFileStream.read(reinterpret_cast<char *>(&value), taglength); 
      offset2 = value;
    } 
    else if(tag==AnglesPhiTag && !m_IsDoppler)
    {
      int len = taglength/sizeof(double);
      double * angles = new double[len];
      inputFileStream.read(reinterpret_cast<char *>(angles), sizeof( double ) * len); 
      double amin = angles[0];
      double amax = angles[len-1];
      double bCentre = (amin+amax)/2;
      for(int i=0; i<len;i++){
        double angle = angles[i]-bCentre;
        this->m_TableAnglesPhi.push_back(std::make_pair(angle, i));
      }
      delete[] angles;
    } 
    else if(tag==AnglesThetaTag && !m_IsDoppler)
    {
      int len = taglength/sizeof(double);
      double * angles = new double[len];
      inputFileStream.read(reinterpret_cast<char *>(angles), sizeof( double ) * len); 
      double amin = angles[0];
      double amax = angles[len-1];
      double bCentre = (amin+amax)/2;
      for(int i=0; i<len;i++){
        double angle = angles[i]-bCentre;
        this->m_TableAnglesTheta.push_back(std::make_pair(angle, i));
      }
      delete[] angles;
    } 

    else if(tag==DimensionXTagDoppler && m_IsDoppler)
    {
      unsigned short dimension;
      inputFileStream.read(reinterpret_cast<char *>(&dimension), taglength); 
      this->SetDimensions(0,dimension);
    } 
    else if(tag==DimensionYTagDoppler && m_IsDoppler)
    {
      unsigned short dimension;
      inputFileStream.read(reinterpret_cast<char *>(&dimension), taglength); 
      this->SetDimensions(1,dimension);
    } 
    else if(tag==DimensionZTagDoppler && m_IsDoppler)
    {
      unsigned short dimension;
      inputFileStream.read(reinterpret_cast<char *>(&dimension), taglength); 
      this->SetDimensions(2,dimension);
    } 
    else if(tag==ResolutionTagDoppler && m_IsDoppler)
    {
      double value;
      inputFileStream.read(reinterpret_cast<char *>(&value), taglength); 
      this->m_Resolution = value *1000.0;
    } 
    else if(tag==Offset1TagDoppler && m_IsDoppler)
    {
      double value;
      inputFileStream.read(reinterpret_cast<char *>(&value), taglength); 
      offset1 = value;
    } 
    else if(tag==Offset2TagDoppler && m_IsDoppler)
    {
      double value;
      inputFileStream.read(reinterpret_cast<char *>(&value), taglength); 
      offset2 = value;
    } 
    else if(tag==AnglesPhiTagDoppler && m_IsDoppler)
    {
      int len = taglength/sizeof(double);
      double * angles = new double[len];
      inputFileStream.read(reinterpret_cast<char *>(angles), sizeof( double ) * len); 
      double amin = angles[0];
      double amax = angles[len-1];
      double bCentre = (amin+amax)/2;
      for(int i=0; i<len;i++){
        double angle = angles[i]-bCentre;
        this->m_TableAnglesPhi.push_back(std::make_pair(angle, i));
      }
      delete[] angles;
    } 
    else if(tag==AnglesThetaTagDoppler && m_IsDoppler)
    {
      int len = taglength/sizeof(double);
      double * angles = new double[len];
      inputFileStream.read(reinterpret_cast<char *>(angles), sizeof( double ) * len); 
      double amin = angles[0];
      double amax = angles[len-1];
      double bCentre = (amin+amax)/2;
      for(int i=0; i<len;i++){
        double angle = angles[i]-bCentre;
        this->m_TableAnglesTheta.push_back(std::make_pair(angle, i));
      }
      delete[] angles;
    } 
    else 
    {
      inputFileStream.seekg(taglength,std::ios::cur);
    }

  }
  m_RadiusBStart = offset1*m_Resolution;
  m_RadiusD = -offset2*m_Resolution;

  inputFileStream.close();

}

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

void KretzImageIO::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  std::string value;
  const MetaDataDictionary & dico = this->GetMetaDataDictionary();
  ExposeMetaData< std::string >(dico, PatientTag, value);
  os << indent << "PatientName: " << value << std::endl;

}
} // end namespace itk
