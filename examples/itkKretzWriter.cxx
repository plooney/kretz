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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkMetaDataObject.h"
#include "itkKretzImageIO.h"
#include "itkKretzImageIOFactory.h"
#include "itkToroidalToCartesianTransform.h"
#include <itkCastImageFilter.h>
#include <itkNormalizeImageFilter.h>
#include <itkResampleImageFilter.h>
#include "itkNearestNeighborInterpolateImageFunction.h"
#include <list>
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int execute(std::string filename, std::string filename_out, std::string filename_cartesian, bool flagNormalise)
{
  typedef unsigned char InputPixelType;
  const unsigned int   Dimension = 3;
  typedef itk::Image< InputPixelType, Dimension > InputImageType;
  typedef itk::Image< double, Dimension > DoubleImageType;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename );
  typedef itk::KretzImageIO           ImageIOType;
  ImageIOType::Pointer kretzImageIO = ImageIOType::New();
  reader->SetImageIO( kretzImageIO );
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "exception in file reader " << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    return EXIT_FAILURE;
  }
  InputImageType::Pointer inputImage = reader->GetOutput();
  //typedef itk::MetaDataDictionary   DictionaryType;
  //DictionaryType & dictionary = inputImage->GetMetaDataDictionary();

  //If a cartesian file is supplied transform it back to toroidal 
  //space using the geometry of the supplied kretz file
  if(!filename_cartesian.empty())
  {
    typedef itk::ImageFileReader< DoubleImageType > ReaderType;
    ReaderType::Pointer reader_cartesian = ReaderType::New();
    reader_cartesian->SetFileName( filename_cartesian );
    reader_cartesian->Update();
    std::cout << filename_cartesian << std::endl;
    typedef itk::ToroidalToCartesianTransform<double,Dimension> T2CTransformType;
    T2CTransformType::Pointer t2c = T2CTransformType::New();
    t2c->SetBModeRadius(kretzImageIO->GetRadiusD());
    t2c->SetSweepRadius(kretzImageIO->GetRadiusBStart());
    t2c->SetResolution(kretzImageIO->GetResolution());
    t2c->SetTablePhi(kretzImageIO->m_TableAnglesPhi);
    t2c->SetTableTheta(kretzImageIO->m_TableAnglesTheta);

    typedef itk::ResampleImageFilter<DoubleImageType,DoubleImageType> ResampleFilterType;
    typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();

    resampleFilter->SetTransform(t2c);
    resampleFilter->SetOutputOrigin(inputImage->GetOrigin());
    resampleFilter->SetOutputDirection(inputImage->GetDirection());
    resampleFilter->SetOutputSpacing(inputImage->GetSpacing());

    typedef itk::NearestNeighborInterpolateImageFunction<DoubleImageType, double >  NearestNeighbourInterpolatorType;
    auto interpolator = NearestNeighbourInterpolatorType::New();
    resampleFilter->SetInterpolator(interpolator);
    resampleFilter->SetSize(inputImage->GetLargestPossibleRegion().GetSize());

    resampleFilter->SetInput(reader_cartesian->GetOutput());
    resampleFilter->Update();
    typedef itk::ImageFileWriter< DoubleImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( resampleFilter->GetOutput() );
    writer->SetFileName( filename_out );
    try
    {
      writer->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      std::cerr << "exception in file writer " << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
      return EXIT_FAILURE;
    }
  } 
  else if(flagNormalise){ //if normalised images wanted outpu double image type

    typedef itk::ImageFileWriter<DoubleImageType> ImageWriterType;
    typedef itk::CastImageFilter<InputImageType,DoubleImageType> CastImageFilterType;
    typedef itk::NormalizeImageFilter<DoubleImageType, DoubleImageType> NormaliseImageFilterType;

    CastImageFilterType::Pointer castImgFilter = CastImageFilterType::New();
    castImgFilter->SetInput(inputImage);
    castImgFilter->Update();

    DoubleImageType::Pointer castImage = castImgFilter->GetOutput();

    NormaliseImageFilterType::Pointer normaliseFilter = NormaliseImageFilterType::New();
    normaliseFilter->SetInput(castImage);
    normaliseFilter->Update();
    DoubleImageType::Pointer normImage = normaliseFilter->GetOutput();

    ImageWriterType::Pointer ITKImageWriter = ImageWriterType::New();

    ITKImageWriter->SetFileName(filename_out.c_str());
    ITKImageWriter->SetInput(normImage);
    ITKImageWriter->Write();

  } 
  else //otherwise write the toroidal volume to the output file 
  {
    typedef itk::ImageFileWriter< InputImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( reader->GetOutput() );
    writer->SetFileName( filename_out );
    try
    {
      writer->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      std::cerr << "exception in file writer " << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{

  std::string filename, filename_cartesian, filename_out;
  po::options_description desc("Allowed options");

  //bool flagMask = false;
  bool flagNormalise = false;
  //bool flagThreshold = false;

  desc.add_options()
    ("help,h", "produce help message")
    ("inputfile,i", po::value<std::string>(&filename)->required(), "set input file")
    ("filename_cartesian,c", po::value<std::string>(&filename_cartesian), "set cartesian file")
    ("outputfile,o", po::value<std::string>(&filename_out)->required(), "set output file")
    ("normalise,n", po::bool_switch(&flagNormalise), "normalise volume")
    ;

  try
  {

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help") || vm.empty() || vm.count("i") || vm.count("o") ) {
      std::cout << desc << "\n";
      return 1;
    }

    po::notify(vm);

    return execute(filename, filename_out, filename_cartesian, flagNormalise);

  }
  catch(std::exception& e)
  {

    std::cerr << "Error: " << e.what() << "\n";

    std::cout << desc << "\n";

    return 1;
  }
  catch(...)
  {
    std::cerr << "Unknown error!" << "\n";
    return 1;
  }

}
