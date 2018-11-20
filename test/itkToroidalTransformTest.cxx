#include <iostream>

#include "itkToroidalToCartesianTransform.h"
#include "itkCartesianToToroidalTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include <itkResampleImageFilter.h>
#include <itkImageFileWriter.h>
#include <cmath>

#include "itkImage.h"
#include "vnl/vnl_math.h"

#include <itkLabelOverlapMeasuresImageFilter.h>

const unsigned int ParametricDimension = 3;

typedef itk::Image<unsigned char, ParametricDimension> ImageType;
typedef itk::LabelOverlapMeasuresImageFilter<ImageType> LabelOverlapType;

std::vector<double> angles1 = {-0.6166550000000001, -0.6097650000000001, -0.602875, -0.5959850000000001, -0.589095, -0.5822050000000001, -0.575315, -0.5684250000000002, -0.5615350000000001, -0.554645, -0.547755, -0.5408650000000002, -0.5339750000000001, -0.527085, -0.520195, -0.5133050000000001, -0.5064150000000001, -0.499525, -0.49263500000000016, -0.4857450000000001, -0.47885500000000003, -0.4719650000000002, -0.4650750000000001, -0.45818500000000006, -0.451295, -0.44440500000000016, -0.4375150000000001, -0.43062500000000004, -0.423735, -0.41684500000000013, -0.40995500000000007, -0.403065, -0.39617500000000017, -0.3892850000000001, -0.38239500000000004, -0.3755050000000002, -0.36861500000000014, -0.3617250000000001, -0.354835, -0.34794499999999995, -0.3410550000000001, -0.33416500000000005, -0.327275, -0.32038500000000014, -0.3134950000000001, -0.306605, -0.2997150000000002, -0.2928250000000001, -0.28593500000000005, -0.2790450000000002, -0.27215500000000015, -0.2652650000000001, -0.258375, -0.25148499999999996, -0.24459500000000012, -0.23770500000000006, -0.230815, -0.22392500000000015, -0.2170350000000001, -0.21014500000000003, -0.20325500000000019, -0.19636500000000012, -0.18947500000000006, -0.18258500000000022, -0.17569500000000016, -0.1688050000000001, -0.16191500000000003, -0.15502499999999997, -0.14813500000000013, -0.14124500000000006, -0.134355, -0.12746500000000016, -0.1205750000000001, -0.11368500000000004, -0.1067950000000002, -0.09990500000000013, -0.09301500000000007, -0.08612500000000023, -0.07923499999999994, -0.0723450000000001, -0.06545500000000026, -0.05856499999999998, -0.05167500000000014, -0.044785000000000075, -0.03789500000000001, -0.03100500000000017, -0.02411500000000011, -0.017225000000000046, -0.010335000000000205, -0.0034450000000001424, 0.0034449999999999203, 0.010334999999999761, 0.017225000000000046, 0.024114999999999887, 0.031004999999999727, 0.03789500000000001, 0.04478499999999985, 0.051674999999999915, 0.05856499999999998, 0.06545499999999982, 0.07234499999999988, 0.07923499999999994, 0.08612499999999979, 0.09301499999999985, 0.09990499999999991, 0.10679499999999975, 0.11368500000000004, 0.12057499999999988, 0.12746499999999972, 0.134355, 0.14124499999999984, 0.1481349999999999, 0.15502499999999997, 0.1619149999999998, 0.16880499999999987, 0.17569499999999993, 0.18258499999999978, 0.18947499999999984, 0.1963649999999999, 0.20325499999999974, 0.21014500000000003, 0.21703499999999987, 0.2239249999999997, 0.230815, 0.23770499999999983, 0.2445949999999999, 0.25148499999999996, 0.2583749999999998, 0.26526499999999986, 0.2721549999999999, 0.27904499999999977, 0.28593499999999983, 0.2928249999999999, 0.29971499999999973, 0.306605, 0.31349499999999986, 0.3203849999999997, 0.327275, 0.3341649999999998, 0.3410549999999999, 0.34794499999999995, 0.3548349999999998, 0.36172499999999985, 0.3686149999999999, 0.37550499999999976, 0.3823949999999998, 0.3892849999999999, 0.3961749999999997, 0.403065, 0.40995499999999985, 0.4168449999999997, 0.423735, 0.4306249999999998, 0.43751499999999965, 0.44440499999999994, 0.4512949999999998, 0.45818500000000006, 0.4650749999999999, 0.47196499999999975, 0.47885500000000003, 0.48574499999999987, 0.4926349999999997, 0.499525, 0.5064149999999998, 0.5133049999999997, 0.520195, 0.5270849999999998, 0.5339749999999996, 0.5408649999999999, 0.5477549999999998, 0.554645, 0.5615349999999999, 0.5684249999999997, 0.575315, 0.5822049999999999, 0.5890949999999997, 0.595985, 0.6028749999999998, 0.6097649999999997, 0.616655};
std::vector<double> angles2 = {-0.7376595572958549, -0.7271963011639988, -0.7167330450321426, -0.7062697889002866, -0.6958065327684305, -0.6853432766365744, -0.6748800205047183, -0.6644167643728622, -0.6539535082410061, -0.64349025210915, -0.6330269959772938, -0.6225637398454378, -0.6121004837135817, -0.6016372275817256, -0.5911739714498695, -0.5807107153180133, -0.5702474591861573, -0.5597842030543012, -0.5493209469224452, -0.5388576907905891, -0.5283944346587328, -0.5179311785268768, -0.5074679223950207, -0.49700466626316464, -0.4865414101313086, -0.4760781539994525, -0.46561489786759624, -0.4551516417357402, -0.44468838560388413, -0.43422512947202807, -0.423761873340172, -0.41329861720831573, -0.4028353610764597, -0.3923721049446036, -0.38190884881274756, -0.3714455926808915, -0.3609823365490352, -0.35051908041717916, -0.3400558242853231, -0.32959256815346705, -0.319129312021611, -0.30866605588975493, -0.2982027997578989, -0.2877395436260426, -0.27727628749418654, -0.2668130313623305, -0.2563497752304744, -0.24588651909861836, -0.23542326296676208, -0.22496000683490602, -0.21449675070304997, -0.2040334945711939, -0.19357023843933785, -0.18310698230748157, -0.17264372617562573, -0.16218047004376945, -0.1517172139119134, -0.14125395778005734, -0.13079070164820106, -0.12032744551634522, -0.10986418938448894, -0.09940093325263288, -0.08893767712077683, -0.07847442098892055, -0.06801116485706471, -0.05754790872520843, -0.04708465259335237, -0.036621396461496314, -0.026158140329640256, -0.015694884197784198, -0.005231628065927918, 0.00523162806592814, 0.015694884197784198, 0.026158140329640256, 0.036621396461496314, 0.047084652593352594, 0.05754790872520843, 0.06801116485706471, 0.07847442098892077, 0.08893767712077683, 0.0994009332526331, 0.10986418938448894, 0.12032744551634522, 0.13079070164820128, 0.14125395778005734, 0.15171721391191362, 0.16218047004376945, 0.17264372617562573, 0.1831069823074818, 0.19357023843933785, 0.2040334945711939, 0.21449675070304997, 0.22496000683490625, 0.23542326296676208, 0.24588651909861836, 0.2563497752304744, 0.2668130313623305, 0.27727628749418654, 0.2877395436260428, 0.2982027997578989, 0.30866605588975493, 0.319129312021611, 0.32959256815346705, 0.3400558242853233, 0.3505190804171794, 0.36098233654903544, 0.3714455926808915, 0.38190884881274756, 0.3923721049446036, 0.4028353610764599, 0.41329861720831595, 0.423761873340172, 0.4342251294720283, 0.44468838560388413, 0.45515164173573996, 0.4656148978675967, 0.4760781539994525, 0.48654141013130836, 0.49700466626316464, 0.5074679223950209, 0.5179311785268768, 0.528394434658733, 0.5388576907905893, 0.5493209469224452, 0.559784203054301, 0.5702474591861573, 0.5807107153180135, 0.5911739714498694, 0.6016372275817257, 0.6121004837135819, 0.6225637398454378, 0.6330269959772941, 0.6434902521091503, 0.6539535082410062, 0.664416764372862, 0.6748800205047183, 0.6853432766365746, 0.6958065327684304, 0.7062697889002867, 0.716733045032143, 0.7271963011639988, 0.7376595572958546};

std::vector<std::pair<double,double>> createTable(std::vector<double> angles){

    std::vector<std::pair<double,double>> table;

    for(unsigned int i=0; i<angles.size(); i++){
        table.push_back(std::make_pair(angles[i], i));
    }

    return table;

}

ImageType::Pointer createImage(ImageType::SizeType size)
{

    ImageType::IndexType start;
    start.Fill(0);

    ImageType::RegionType region(start, size);

    ImageType::SpacingType spacing;
    spacing.Fill(1);

    ImageType::PointType origin;
    origin.Fill(0);

    ImageType::DirectionType direction;
    direction.SetIdentity();

    ImageType::Pointer returnImage = ImageType::New();
    returnImage->SetRegions(region);
    returnImage->SetOrigin(origin);
    returnImage->SetDirection(direction);
    returnImage->SetSpacing(spacing);
    returnImage->Allocate();
    returnImage->FillBuffer(1);

    for(unsigned int r = 10; r < size[0]-10; r++)
    {
        for(unsigned int c = 10; c < size[1]-10; c++)
        {
            for(unsigned int s = 10; s < size[2]-10; s++)
            {
                ImageType::IndexType pixelIndex;
                pixelIndex[0] = r;
                pixelIndex[1] = c;
                pixelIndex[2] = s;

                returnImage->SetPixel(pixelIndex, 2);
            }
        }
    }


    for(unsigned int r = 40; r < 50; r++)
    {
        for(unsigned int c = 40; c < 50; c++)
        {
            for(unsigned int s = 40; s < 50; s++)
            {
                ImageType::IndexType pixelIndex;
                pixelIndex[0] = r;
                pixelIndex[1] = c;
                pixelIndex[2] = s;

                returnImage->SetPixel(pixelIndex, 3);
            }
        }
    }

    for(unsigned int r = size[0]-50; r < size[0]-40; r++)
    {
        for(unsigned int c = size[1]-50; c <size[1]-40; c++)
        {
            for(unsigned int s = size[2]-50; s < size[2]-40; s++)
            {
                ImageType::IndexType pixelIndex;
                pixelIndex[0] = r;
                pixelIndex[1] = c;
                pixelIndex[2] = s;

                returnImage->SetPixel(pixelIndex, 4);
            }
        }
    }

    return returnImage;
}

int main(int ,char *[] )
{
    const unsigned int Dimension = 3;

    //Test is itkToroidalToCartesian after itkCartesianToToroidal identity

    std::vector<std::pair<double, double> > TableAngles1 = createTable(angles1);
    std::vector<std::pair<double, double> > TableAngles2 = createTable(angles2);

    typedef itk::ToroidalToCartesianTransform<double,Dimension> T2CTransformType;
    typedef itk::CartesianToToroidalTransform<double,Dimension> C2TTransformType;
    typedef itk::ImageFileWriter<ImageType> ImageWriterType;

    T2CTransformType::Pointer t2c = T2CTransformType::New();
    C2TTransformType::Pointer c2t = C2TTransformType::New();

    ImageType::DirectionType direction;
    direction.SetIdentity();

    double rBstart=41.599998474121094;
    double rD=-19.170000076293945;
    double rResol =0.20533333718776703;

    ImageType::SpacingType cartesianSpacing;
    cartesianSpacing.Fill(0.6);

    c2t->SetBModeRadius(rD);
    c2t->SetSweepRadius(rBstart);
    c2t->SetResolution(rResol);
    c2t->SetTableTheta(TableAngles1);
    c2t->SetTablePhi(TableAngles2);

    t2c->SetBModeRadius(rD);
    t2c->SetSweepRadius(rBstart);
    t2c->SetResolution(rResol);
    t2c->SetTableTheta(TableAngles1);
    t2c->SetTablePhi(TableAngles2);

    ImageType::PointType p1;
    p1[0] = -10;
    p1[1] = 10;
    p1[2] = 10;
    ImageType::PointType p2 = c2t->TransformPoint(p1);
    ImageType::PointType p3 = t2c->TransformPoint(p2);
    std::cout << p1 << p2 << p3 << std::endl;

    ImageType::PointType cartesianOrigin;
    cartesianOrigin[0] = -122.81173706054688;
    cartesianOrigin[1] = -125.04386901855469;
    cartesianOrigin[2] = 20.098983764648438;

    ImageType::SizeType cartesian_size = {410, 418, 318};
    ImageType::SizeType size_toroidal = {580, 180, 142};

    typedef itk::ResampleImageFilter<ImageType,ImageType> ResampleFilterType;
    typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();

    resampleFilter->SetTransform(c2t);
    resampleFilter->SetOutputOrigin(cartesianOrigin);
    resampleFilter->SetOutputDirection(direction);
    resampleFilter->SetOutputSpacing(cartesianSpacing);

    typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double >  NearestNeighbourInterpolatorType;
    auto interpolator = NearestNeighbourInterpolatorType::New();
    resampleFilter->SetInterpolator(interpolator);
    resampleFilter->SetSize(cartesian_size);

    ImageType::Pointer toroidalImage = createImage( size_toroidal);
    resampleFilter->SetInput(toroidalImage);
    resampleFilter->Update();
    ImageType::Pointer cartesianImage = resampleFilter->GetOutput();

    ImageWriterType::Pointer ITKImageWriter = ImageWriterType::New();
    ITKImageWriter->SetFileName("test/testImageCartesian.nii.gz");
    ITKImageWriter->SetInput(cartesianImage);
    ITKImageWriter->Write();

    resampleFilter = ResampleFilterType::New();
    resampleFilter->SetTransform(t2c);

    ImageType::PointType origin;
    origin.Fill(0);

    ImageType::SizeType size_toroidal_output;
    size_toroidal_output[0] = size_toroidal[0];
    size_toroidal_output[1] = size_toroidal[1];
    size_toroidal_output[2] = size_toroidal[2];

    ImageType::SpacingType spacing;
    spacing.Fill(1);

    resampleFilter->SetOutputOrigin(origin);
    resampleFilter->SetOutputDirection(direction);
    resampleFilter->SetOutputSpacing(spacing);

    typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double >  NearestNeighbourInterpolatorType;
    resampleFilter->SetInterpolator(interpolator);
    resampleFilter->SetSize(size_toroidal);
    resampleFilter->SetInput(cartesianImage);
    resampleFilter->Update();

    ImageType::Pointer toroiidalFromCartesianImage = resampleFilter->GetOutput();
    ITKImageWriter->SetFileName("test/testImageToroidal.nii.gz");
    ITKImageWriter->SetInput(toroiidalFromCartesianImage);
    ITKImageWriter->Write();

    resampleFilter = ResampleFilterType::New();
    resampleFilter->SetTransform(c2t);
    resampleFilter->SetOutputOrigin(cartesianOrigin);
    resampleFilter->SetOutputDirection(direction);
    resampleFilter->SetOutputSpacing(cartesianSpacing);
    resampleFilter->SetInterpolator(interpolator);
    resampleFilter->SetSize(cartesian_size);
    resampleFilter->SetInput(toroiidalFromCartesianImage);
    resampleFilter->Update();

    ImageType::Pointer cartesianImage2 = resampleFilter->GetOutput();
    ITKImageWriter->SetFileName("test/testImageCartesian2.nii.gz");
    ITKImageWriter->SetInput(cartesianImage2);
    ITKImageWriter->Write();

    LabelOverlapType::Pointer metric = LabelOverlapType::New();
    metric->SetSourceImage(cartesianImage);
    metric->SetTargetImage(cartesianImage2);
    metric->Update();

    double fn_error = metric->GetFalseNegativeError();
    double fp_error = metric->GetFalsePositiveError();
    double overlap1 = metric->GetTargetOverlap(1);
    double overlap2 = metric->GetTargetOverlap(2);
    double overlap3 = metric->GetTargetOverlap(3);
    double overlap4 = metric->GetTargetOverlap(4);
    double mean_overlap = metric->GetMeanOverlap();

    std::cout << "FN error: " << fn_error << std::endl;
    std::cout << "FP error: " << fp_error << std::endl;
    std::cout << "Overlap 1: " << overlap1 << std::endl;
    std::cout << "Overlap 2: " << overlap2 << std::endl;
    std::cout << "Overlap 3: " << overlap3 << std::endl;
    std::cout << "Overlap 4: " << overlap4 << std::endl;
    std::cout << "Mean overlap: " << mean_overlap << std::endl;

    if(mean_overlap < 0.99) {
        std::cout << "Failed overlap test" << std::endl;
        return EXIT_FAILURE;
    } else {
        std::cout << "Passed overlap test" << std::endl;
        return EXIT_SUCCESS;
    }


}
