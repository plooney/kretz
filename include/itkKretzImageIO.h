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
};
static Tag PatientTag(0x0110, 0x0002);
static Tag DimensionXTag(0xc000, 0x0001);
static Tag DimensionYTag(0xc000, 0x0002);
static Tag DimensionZTag(0xc000, 0x0003);
static Tag ResolutionTag(0xc100, 0x0001);
static Tag Offset1Tag(0xc200, 0x0001);
static Tag Offset2Tag(0xc200, 0x0002);
static Tag Angles1Tag(0xc300, 0x0001);
static Tag Angles2Tag(0xc300, 0x0002);
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
 * \ingroup ITKIOGDCM
 *
 * \wiki
 * \wikiexample{DICOM/ResampleDICOM,Resample a DICOM series}
 * \endwiki
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

  /** Set/Get the original component type of the image. This differs from
   * ComponentType which may change as a function of rescale slope and
   * intercept. */
  itkGetEnumMacro(InternalComponentType, IOComponentType);
  itkSetEnumMacro(InternalComponentType, IOComponentType);

  /*-------- This part of the interfaces deals with writing data. ----- */

  /** Determine the file type. Returns true if this ImageIO can write the
   * file specified. GDCM triggers on ".dcm" and ".dicom". */
  virtual bool CanWriteFile(const char *) ITK_OVERRIDE;

  /** Writes the spacing and dimensions of the image.
   * Assumes SetFileName has been called with a valid file name. */
  virtual void WriteImageInformation() ITK_OVERRIDE;

  /** Writes the data to disk from the memory buffer provided. Make sure
   * that the IORegion has been set properly. */
  virtual void Write(const void *buffer) ITK_OVERRIDE;

  void GetPatientName(char *name);

  /** Set/Get a compression type to use. */
  typedef enum { JPEG = 0, JPEG2000, JPEGLS, RLE } TCompressionType;
  itkSetEnumMacro(CompressionType, TCompressionType);
  itkGetEnumMacro(CompressionType, TCompressionType);

  typedef std::vector<std::pair<double, double>> TTableAngleType;
  std::vector<std::pair<double, double>> m_TableAngles1;
  std::vector<std::pair<double, double>> m_TableAngles2;

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

  double m_Offset1;
  double m_Offset2;
  double m_Resolution;
  double m_rBstart;
  double m_rD;



private:
  ITK_DISALLOW_COPY_AND_ASSIGN(KretzImageIO);

  std::string m_PatientName;

  /** defines whether this image is a 2D out of a 2D image
   *  or a 2D out of a 3D image. */
  unsigned int     m_GlobalNumberOfDimensions;
  TCompressionType m_CompressionType;

  ImageIOBase::IOComponentType m_InternalComponentType;
};
} // end namespace itk

#endif //itkKretzImageIO_h
