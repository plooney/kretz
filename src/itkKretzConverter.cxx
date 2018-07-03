#include "stdlib.h"
#include <sstream>
#include <iostream>
#include <string>

#include "itkImage.h"

#include "itkResampleImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkCartesianToToroidalTransform.h"
#include "itkToroidalToCartesianTransform.h"
#include "itkKretzImageIO.h"
#include <itkCastImageFilter.h>
#include <itkNormalizeImageFilter.h>
#include <itkInterpolateImageFunction.h>

#include <itkPointSetToImageFilter.h>
#include <itkMeshSpatialObject.h>
#include <itkBoundingBox.h>

#include <ctime>
#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <itkRescaleIntensityImageFilter.h>

#include <unistd.h>    /* for getopt */

using namespace std;
namespace po = boost::program_options;

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> ImageType;
typedef itk::Image<double, Dimension> DoubleImageType;
typedef itk::ResampleImageFilter<ImageType,ImageType> ResampleImageFilterType;
typedef itk::ImageFileWriter<ImageType> ImageWriterType;
typedef itk::LinearInterpolateImageFunction<ImageType, double> LinearInterpolatorType;
typedef itk::IdentityTransform<double, Dimension> TransformType;
typedef itk::ImageFileReader<ImageType> ImageFileReaderType;
typedef itk::ToroidalToCartesianTransform<double, Dimension> T2CTransformType;
typedef itk::CartesianToToroidalTransform<double, Dimension> C2TTransformType;
typedef itk::Mesh< double, 3 > MeshType;
typedef itk::PointSet<ImageType::PixelType, Dimension> PointSetType;
typedef PointSetType::PointsContainerPointer PointsContainerPointer;
typedef itk::BoundingBox<PointSetType::PointIdentifier, Dimension, double, PointSetType::PointsContainer> BoundingBoxType;

/*
 * Used to find a binary mask of a toroidal volume in cartesian coordinates
 */
ImageType::Pointer createMaskImage(ImageType::Pointer image)
{
  ImageType::Pointer returnImage = ImageType::New();
  returnImage->SetOrigin(image->GetOrigin());
  returnImage->SetDirection(image->GetDirection());
  returnImage->SetSpacing(image->GetSpacing());
  returnImage->SetRegions(image->GetLargestPossibleRegion());
  returnImage->Allocate();
  returnImage->FillBuffer(1);

  return returnImage;
}

/*
 * We need to find the bounds of the toroidal volume in cartesian coordinates 
 * We cannot assume the transformed toroidal boundary points form the cartesian bounds 
 * Hence we need to transform each point and find the bounds of the computed mesh
 */
BoundingBoxType::BoundsArrayType computeBounds(ImageType::Pointer image, T2CTransformType::Pointer t2c)
{
  MeshType::Pointer mesh = MeshType::New();
  itk::ImageRegionConstIteratorWithIndex<ImageType> imageIterator(image,image->GetLargestPossibleRegion());
  while(!imageIterator.IsAtEnd()) 
  {
    ImageType::IndexType index = imageIterator.GetIndex();
    ImageType::PointType p;
    p[0] = index[0];
    p[1] = index[1];
    p[2] = index[2];
    ImageType::PointType p_cart = t2c->TransformPoint(p); 
    
    unsigned long i = mesh->GetNumberOfPoints();
    mesh->SetPoint(i, p_cart);
    ++imageIterator;
  }
  PointsContainerPointer points = mesh->GetPoints();

  BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();

  boundingBox->SetPoints(points);
  boundingBox->ComputeBoundingBox();

  BoundingBoxType::BoundsArrayType bounds = boundingBox->GetBounds();

  return bounds;
}


int execute(std::string filename, std::string filename_out, std::vector<int> size_vec,std::vector<float> resol_vec, bool flagMask, bool flagNormalise, bool flagThreshold, bool isDoppler)
{
  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( filename );

  typedef itk::KretzImageIO           ImageIOType;
  ImageIOType::Pointer kretzImageIO = ImageIOType::New();
  reader->SetImageIO( kretzImageIO );
  reader->Update();

  ImageType::Pointer toroidalImage = reader->GetOutput();
  T2CTransformType::Pointer t2c = T2CTransformType::New();
  t2c->SetBModeRadius(kretzImageIO->GetrD());
  t2c->SetSweepRadius(kretzImageIO->GetrBstart());
  t2c->SetResolution(kretzImageIO->GetResolution());
  t2c->SetTableTheta(kretzImageIO->m_TableAnglesTheta);
  t2c->SetTablePhi(kretzImageIO->m_TableAnglesPhi);
  

  //Find the bounds of the toroidal volume in cartesian coordinates
  BoundingBoxType::BoundsArrayType bounds = computeBounds(toroidalImage, t2c);
  if(isDoppler){
	  reader = ReaderType::New();
	  reader->SetFileName( filename );
	  kretzImageIO = ImageIOType::New();
	  kretzImageIO->SetisDoppler(isDoppler);
	  reader->SetImageIO( kretzImageIO );
	  reader->Update();
	  toroidalImage = reader->GetOutput();
  }

  std::cout << "size " << size_vec.size() << std::endl;
  std::cout << "resol " << resol_vec.size() << std::endl;

  if(size_vec[0]==0) //if the resolution is provided set the size
  {
	  size_vec[0]=(int) ((bounds[1]-bounds[0])/resol_vec[0]);
	  size_vec[1]=(int) ((bounds[3]-bounds[2])/resol_vec[1]);
	  size_vec[2]=(int) ((bounds[5]-bounds[4])/resol_vec[2]);
  } 
  else if(resol_vec[0]==0) //if the size is provided set the resolution
  {
	  resol_vec[0]=(bounds[1]-bounds[0])/(size_vec[0]-1);
	  resol_vec[1]=(bounds[3]-bounds[2])/(size_vec[1]-1);
	  resol_vec[2]=(bounds[5]-bounds[4])/(size_vec[2]-1);
  }
  else
  {
    std::cerr << "Error: both resolution and size provided" <<  std::endl;

    return EXIT_FAILURE;
  }
  std::cout << "resol " << resol_vec.at(0) << " " << resol_vec.at(1) << " " << resol_vec.at(2) << std::endl;
  std::cout << "size " << size_vec[0] << " " << size_vec[1] << " " << size_vec[2]  << std::endl;

  //if a mask is required create a binary image with ones everywhere 
  //and same dimensions as input file
  if(flagMask) toroidalImage = createMaskImage(toroidalImage);

  if(flagNormalise){ //if normalised images wanted outpu double image type

	  typedef itk::ImageFileWriter<DoubleImageType> ImageWriterType;
	  typedef itk::CastImageFilter<ImageType,DoubleImageType> CastImageFilterType;
	  typedef itk::NormalizeImageFilter<DoubleImageType, DoubleImageType> NormaliseImageFilterType;

	  CastImageFilterType::Pointer castImgFilter = CastImageFilterType::New();
	  castImgFilter->SetInput(toroidalImage);
	  castImgFilter->Update();

	  DoubleImageType::Pointer castImage = castImgFilter->GetOutput();

	  NormaliseImageFilterType::Pointer normaliseFilter = NormaliseImageFilterType::New();
	  normaliseFilter->SetInput(castImage);
	  normaliseFilter->Update();
	  DoubleImageType::Pointer normImage = normaliseFilter->GetOutput();

	  DoubleImageType::PointType origin;
	  origin[0] = bounds[0];
	  origin[1] = bounds[2];
	  origin[2] = bounds[4];

	  C2TTransformType::Pointer c2t = C2TTransformType::New();
	  c2t->SetBModeRadius(kretzImageIO->GetrD());
	  c2t->SetSweepRadius(kretzImageIO->GetrBstart());
	  c2t->SetResolution(kretzImageIO->GetResolution());
	  c2t->SetTableTheta(kretzImageIO->m_TableAnglesTheta);
	  c2t->SetTablePhi(kretzImageIO->m_TableAnglesPhi);
	  typedef itk::ResampleImageFilter<DoubleImageType,DoubleImageType> ResampleFilterType;

	  typedef DoubleImageType::SpacingType SpacingType;
	  SpacingType spacing;
	  spacing[0] = resol_vec.at(0);
	  spacing[1] = resol_vec.at(1);
	  spacing[2] = resol_vec.at(2);

	  typedef typename DoubleImageType::SizeType SizeType;
	  SizeType size;
	  size[0]= size_vec[0];
	  size[1]= size_vec[1];
	  size[2]= size_vec[2];

	  typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
          resampleFilter->SetInput(normImage);
	  resampleFilter->SetTransform(c2t);
	  resampleFilter->SetSize(size);
	  resampleFilter->SetOutputOrigin(origin);
	  resampleFilter->SetOutputSpacing(spacing);
	  resampleFilter->Update();
	  DoubleImageType::Pointer output = resampleFilter->GetOutput(); 

	  ImageWriterType::Pointer ITKImageWriter = ImageWriterType::New();

	  ITKImageWriter->SetFileName(filename_out.c_str());
	  ITKImageWriter->SetInput(output);
	  ITKImageWriter->Write();

  } 
  else //If non normalised images required can revert to original image type
  {
	  ImageType::PointType origin;
	  origin[0] = bounds[0];
	  origin[1] = bounds[2];
	  origin[2] = bounds[4];

	  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
	  C2TTransformType::Pointer c2t = C2TTransformType::New();
	  c2t->SetBModeRadius(kretzImageIO->GetrD());
	  c2t->SetSweepRadius(kretzImageIO->GetrBstart());
	  c2t->SetResolution(kretzImageIO->GetResolution());
	  c2t->SetTableTheta(kretzImageIO->m_TableAnglesTheta);
	  c2t->SetTablePhi(kretzImageIO->m_TableAnglesPhi);
	  typedef itk::ResampleImageFilter<ImageType,ImageType> ResampleFilterType;

	  typedef ImageType::SpacingType SpacingType;
	  SpacingType spacing;
	  spacing[0] = resol_vec[0];
	  spacing[1] = resol_vec[1];
	  spacing[2] = resol_vec[2];

	  typedef typename ImageType::SizeType SizeType;
	  SizeType size;
	  size[0]= size_vec[0];
	  size[1]= size_vec[1];
	  size[2]= size_vec[2];


	  typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
    resampleFilter->SetInput(toroidalImage);

	  resampleFilter->SetTransform(c2t);
	  resampleFilter->SetSize(size);
	  resampleFilter->SetOutputOrigin(origin);
	  resampleFilter->SetOutputSpacing(spacing);
	  resampleFilter->Update();
	  ImageType::Pointer output = resampleFilter->GetOutput(); 

	  ImageWriterType::Pointer ITKImageWriter = ImageWriterType::New();

	  ITKImageWriter->SetFileName(filename_out.c_str());
	  ITKImageWriter->SetInput(output);
	  ITKImageWriter->Write();
  }
  return EXIT_SUCCESS;

}

int main(int argc, char ** argv)
{

    std::string filename, filename_out;
    std::vector<int> size_vec(3);
    std::vector<float> resol_vec(3);

    po::options_description desc("Allowed options");

    bool flagMask = false;
    bool flagNormalise = false;
    bool flagThreshold = false;
    bool flagDoppler = false;

    desc.add_options()
            ("help,h", "produce help message")
            ("inputfile,i", po::value<std::string>(&filename)->required(), "set input file")
            ("outputfile,o", po::value<std::string>(&filename_out)->required(), "set output file")
            ("resol,r", po::value<std::vector<float>>(&resol_vec)->multitoken(), "set resolution")
            ("size,s", po::value<std::vector<int>>(&size_vec)->multitoken(), "set size")
            ("mask,m", po::bool_switch(&flagMask), "mask volume")
            ("normalise,n", po::bool_switch(&flagNormalise), "normalise volume")
            ("threshold,t", po::bool_switch(&flagThreshold), "threshold volume")
            ("isDoppler,d", po::bool_switch(&flagDoppler), "output power Doppler")
            ;

    try
    {

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help") || vm.empty() || (vm.count("size")==0 && vm.count("resol")==0)) {
            cout << desc << "\n";
            return EXIT_FAILURE;
        }

        po::notify(vm);

        return execute(filename,filename_out,size_vec,resol_vec,flagMask,flagNormalise,flagThreshold,flagDoppler);

    }
    catch(std::exception& e)
    {

        std::cerr << "Error: " << e.what() << "\n";

        cout << desc << "\n";

        return EXIT_FAILURE;
    }
    catch(...)
    {
        std::cerr << "Unknown error!" << "\n";
        return EXIT_FAILURE;
    }

}
