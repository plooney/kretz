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


/** \class KretzImageIO
 *
 *
 *  \brief ImageIO class for reading Kretzfile V1.0
 *
 * The pixel spacing in the toroidal space is not fixed. The variable angel
 * spacing and radial resolution are stored in class variables. Hence, the
 * spacing and origin do not correspond to real world values.
 *
 * Writing is not supported.
 *
 * \ingroup IOFilters
 *
 * \ingroup IOKretz
 *
 */
class KretzImageIO:public ImageIOBase
{
public:
  /** Standard class typedefs. */
  typedef KretzImageIO         Self;
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

  typedef std::vector<std::pair<double, double> > TTableAngleType;
  std::vector<std::pair<double, double> > m_TableAnglesTheta;
  std::vector<std::pair<double, double> > m_TableAnglesPhi;

  itkSetMacro(RadiusBStart, double);
  itkGetMacro(RadiusBStart, double);
  itkGetMacro(RadiusD, double);
  itkSetMacro(RadiusD, double);
  itkGetMacro(Resolution, double);
  itkSetMacro(Resolution, double);
  itkGetMacro(IsDoppler, bool);
  itkSetMacro(IsDoppler, bool);

protected:
  KretzImageIO();
  ~KretzImageIO() ITK_OVERRIDE;
  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  double m_Resolution;
  double m_RadiusBStart;
  double m_RadiusD;

  bool m_IsDoppler = false;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(KretzImageIO);

  std::string m_PatientName;

};
} // end namespace itk

#endif //itkKretzImageIO_h
