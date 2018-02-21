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
#ifndef itkKretzImageIO_h
#define itkKretzImageIO_h

#include "itkImageIOBase.h"
#include <fstream>
#include <string>

namespace itk
{


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
static Tag DimensionXTag(0xc000, 0x0001);
static Tag DimensionYTag(0xc000, 0x0002);
static Tag DimensionZTag(0xc000, 0x0003);
static Tag ResolutionTag(0xc100, 0x0001);
static Tag Offset1Tag(0xc200, 0x0001);
static Tag Offset2Tag(0xc200, 0x0002);
static Tag AnglesPhiTag(0xc300, 0x0001);
static Tag AnglesThetaTag(0xc300, 0x0002);
static Tag ImageTag(0xd000, 0x0001);
static Tag CineFramesTag(0xd400, 0x0001);
static Tag SizeFramesTag(0xd400, 0x0002);
static Tag TimingFramesTag(0xd400, 0x0005);
static Tag Image4dTag(0xd600, 0x0001);

/** \class KretzImageIO
 *
 *  \brief ImageIO class for reading Kretzfile V1.0 
 *  The pixel spacing in the toroidal space is not fixed. The variable angel spacing and radial resolution are stored in class variables 
 *  Writing is not supported
 *
 *  \ingroup IOFilters
 *
 * \ingroup ITKIOKRETZ
 *
 */
class InternalHeader;
class KretzImageIO:public ImageIOBase
{
public:
  /** Standard class typedefs. */
  typedef KretzImageIO          Self;
  typedef ImageIOBase          Superclass;
  typedef SmartPointer< Self > Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(KretzImageIO, Superclass);

  /*-------- This part of the interface deals with reading data. ------ */

  /** Determine the file type. Returns true if this ImageIO can read the
   * file specified. */
  virtual bool CanReadFile(const char *) ITK_OVERRIDE;

  /** Set the spacing and dimesion information for the current filename. */
  virtual void ReadImageInformation() ITK_OVERRIDE;

  /** Reads the data from disk into the memory buffer provided. */
  virtual void Read(void *buffer) ITK_OVERRIDE;

  /*-------- This part of the interfaces deals with writing data. ----- */

  //Writing not supported
  virtual bool CanWriteFile(const char *) ITK_OVERRIDE;

  //Writing not supported
  virtual void WriteImageInformation() ITK_OVERRIDE;

  //Writing not supported
  virtual void Write(const void *buffer) ITK_OVERRIDE;

  void GetPatientName(char *name);

  typedef std::vector<std::pair<double, double>> TTableAngleType;
  std::vector<std::pair<double, double>> m_TableAnglesTheta;
  std::vector<std::pair<double, double>> m_TableAnglesPhi;

  itkSetMacro(rBstart, double);
  itkGetMacro(rBstart, double);
  itkGetMacro(rD, double);
  itkSetMacro(rD, double);
  itkGetMacro(Resolution, double);
  itkSetMacro(Resolution, double);


protected:
  KretzImageIO();
  ~KretzImageIO() ITK_OVERRIDE;
  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  double m_Resolution;
  double m_rBstart;
  double m_rD;



private:
  ITK_DISALLOW_COPY_AND_ASSIGN(KretzImageIO);

  std::string m_PatientName;

};
} // end namespace itk

#endif //itkKretzImageIO_h
