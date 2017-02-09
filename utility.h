//Template for the Dimension
//Image Properties
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#ifndef utility_h
#define utility_h

#include "itkImage.h"
#include <itkIndex.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include <itkNiftiImageIO.h>

#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>


#include "itkTimeProbe.h"
itk::TimeProbe t;
std:: ofstream myfile2;
const unsigned int  Dimension = 3;
typedef  float InternalPixelType;
typedef short int OutputPixelType;

typedef itk::Image<InternalPixelType, Dimension> InternalImageType;

//typedef unsigned char OutputPixelType;
typedef itk::Image<OutputPixelType, Dimension> OutputImageType;
typedef  InternalImageType::IndexType  IndexType;
typedef  InternalImageType::IndexType  Indexarray[]; 

typedef long unsigned int WatershedPixelType;
typedef itk::Image<WatershedPixelType, Dimension> WatershedImageType;

InternalImageType::Pointer imageSeg;
InternalImageType::RegionType region;
InternalImageType::SizeType Imagesize;

using namespace std;

typedef itk::ImageFileReader<InternalImageType> ReaderType;  //the diacom image
ReaderType::Pointer readImage;
typedef itk::ImageFileReader<OutputImageType> ReaderTypeout;  //the diacom image
ReaderTypeout::Pointer readImageOut;


static InternalImageType::Pointer ReadImage(char *filename)
{
	readImage = ReaderType::New();
	readImage->SetFileName(filename);
	try
	{
		readImage->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e << std::endl;
		return NULL;
	}
	imageSeg = readImage->GetOutput();
	region = imageSeg->GetLargestPossibleRegion();
	Imagesize = region.GetSize();
	myfile2 << Imagesize << ",";
	return imageSeg;
}
static OutputImageType::Pointer ReadImageOut(char *filename)
{
	readImageOut = ReaderTypeout::New();
	readImageOut->SetFileName(filename);
	try
	{
		readImageOut->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e << std::endl;
		return NULL;
	}

	//std::cout << "Input image size: " << readImageOut->GetOutput()->GetBufferedRegion().GetSize() << std::endl;
	return readImageOut->GetOutput();;
}
static int Write(InternalImageType::Pointer image, std::string parameter)
{
	typedef itk::ImageFileWriter<InternalImageType> WriterType;
	WriterType::Pointer writeImage = WriterType::New();

	writeImage->SetInput(image);
	writeImage->SetFileName(parameter);
//	std::cout << "Writting Image:: " << parameter << endl;
	myfile2 << "Writting Image:: " << parameter << endl;
	try
	{
		writeImage->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception caught here!! !" << std::endl;
		std::cerr << e << std::endl;
		return 0;
	}
	return 1;
}
static int Write(OutputImageType::Pointer image, std::string parameter)
{
	typedef itk::ImageFileWriter<OutputImageType> WriterType;
	WriterType::Pointer writeImage = WriterType::New();
	writeImage->SetInput(image);
	writeImage->SetFileName(parameter);
	//cout << "Writting Image:: " << parameter << endl;
	myfile2 << "Writting Image:: " << parameter << endl;
	try
	{
		writeImage->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << e << std::endl;
		return 0;
	}
	return 1;
}
static int Write(WatershedImageType::Pointer image, std::string parameter)
{
	typedef itk::ImageFileWriter<WatershedImageType> WriterType;
	WriterType::Pointer writeImage = WriterType::New();
	writeImage->SetInput(image);
	writeImage->SetFileName(parameter);
	//cout << "Writting Image:: " << parameter << endl;
	myfile2 << "Writting Image:: " << parameter << endl;
	try
	{
		writeImage->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << e << std::endl;
		return 0;
	}
	return 1;
}

static OutputImageType::Pointer normalizeGtImage(InternalImageType::Pointer image)
{
	typedef itk::RescaleIntensityImageFilter<InternalImageType, OutputImageType> RescaleFilterType;
	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	//Image GrayValue Limits
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(1);
	rescaler->SetInput(image);
	//rescaler->SetInPlace(true);
	//std::cout << "normalizing Image between 0-1..." << std::endl;
	//myfile2 << "normalizing Image between 0-1..." << std::endl;
	//Rescale the Image
	t.Start();
	rescaler->Update();
	t.Stop();
	//std::cout << "...normalize Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	//myfile2 << "...normalize Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;

	return rescaler->GetOutput();
}
static InternalImageType::Pointer RescaleImage(InternalImageType::Pointer image)
{
	typedef itk::RescaleIntensityImageFilter<InternalImageType, InternalImageType> RescaleFilterType;
	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	//Image GrayValue Limits
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	//Anisotropic Diffusion Output as Rescaler Input
	rescaler->SetInput(image);
	rescaler->SetInPlace(true);
	//cout << image->GetLargestPossibleRegion().GetSize() << endl;
	//cout << image->GetBufferedRegion().GetSize() << endl;
//	std::cout<<"Rescaling Image between 0-255..."<<std::endl;
	//myfile2<< "Rescaling Image between 0-255..." << std::endl;
	//Rescale the Image
	t.Start();
	rescaler->Update();
	t.Stop();
	//std::cout << "...Rescale Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	//myfile2<< "...Rescale Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	//Write(rescaler->GetOutput(), "rescaled.nii.gz");
	return rescaler->GetOutput();
}
static float* evaluate(OutputImageType::Pointer imageGt, OutputImageType::Pointer imageSeg)
{
	typedef itk::LabelOverlapMeasuresImageFilter<OutputImageType> FilterType;
	FilterType::Pointer labelFilter = FilterType::New();
	labelFilter->SetSourceImage(imageGt);
	labelFilter->SetTargetImage(imageSeg);
	try
	{
		labelFilter->Update();
	}
	catch (itk::ExceptionObject & e) {
		std::cerr << "label Filter Error" << std::endl;
		std::cerr << e << std::endl;
	}
	std::cout << std::setw(17) << "Union (jaccard)" << std::setw(17) << labelFilter->GetJaccardCoefficient() << std::endl;
	std::cout << std::setw(17) << "Mean (dice)" << std::setw(17) << labelFilter->GetDiceCoefficient() << std::endl;
	std::cout << std::setw(17) << "False negative" << std::setw(17) << labelFilter->GetFalseNegativeError() << std::endl;
	std::cout << std::setw(17) << "False positive" << std::setw(17) << labelFilter->GetFalsePositiveError() << std::endl;
	float* vec = 0;
	vec = new float[3];
	vec[0] = labelFilter->GetJaccardCoefficient();
    vec[1] = labelFilter->GetDiceCoefficient();
	vec[2] = labelFilter->GetFalseNegativeError();
	vec[3] = labelFilter->GetFalsePositiveError();
	return vec;
}
static float getStats(InternalImageType::Pointer imageStats)
{
	// Compute some statistics
	typedef itk::StatisticsImageFilter<InternalImageType> StatisticsImageFilterType;
	StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
	statisticsImageFilter->SetInput(imageStats);
	statisticsImageFilter->Update();
	float mean = statisticsImageFilter->GetMean();
	float sigma = statisticsImageFilter->GetSigma();
	// Print statistics
	std::cout << std::endl;
	std::cout << "----------------------------------------------------" << std::endl;
	std::cout << "Image maximum: " << statisticsImageFilter->GetMaximum() << std::endl;
	std::cout << "Image maximum: " << statisticsImageFilter->GetMinimum() << std::endl;
	std::cout << "Mean of image: " << mean << std::endl;
	std::cout << "Sigma of image: " << sigma << std::endl;
	std::cout << "----------------------------------------------------" << std::endl;
	return 0;
}
static float getVolume(OutputImageType::Pointer imageStats, double sourceImageSpacing[])
{
	// Compute some statistics
	typedef itk::StatisticsImageFilter<OutputImageType> StatisticsImageFilterType;
	StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
	statisticsImageFilter->SetInput(imageStats);
	statisticsImageFilter->Update();
	float volume = (statisticsImageFilter->GetSum())*sourceImageSpacing[0] * sourceImageSpacing[1] * sourceImageSpacing[2];
	/*cout << volume << " mm³" << std::endl;
	std::cout << "----------------------------------------------------" << std::endl;
	myfile2 << volume << " mm³" << std::endl;
	myfile2 << "----------------------------------------------------" << std::endl;*/
	return volume;
}
static float* getMuSigma(InternalImageType::Pointer image)
{
	float intensity = 0.0;
	int numInteriorPixels1 = 0;
	itk::ImageRegionIterator<InternalImageType>  it(image, image->GetLargestPossibleRegion());
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		if (image->GetPixel(it.GetIndex()))
		{
			numInteriorPixels1++;
			intensity += image->GetPixel(it.GetIndex());
		}
	}
	float mean = 0.0, sigma = 0.0;
	mean = intensity / numInteriorPixels1;
	// std::cout << "Mean of the region is "<<mean <<" "<< std::endl;
	double variance = 0;
	// numInteriorPixels1 = 0;
	for (it.GoToBegin(); !it.IsAtEnd(); ++it) // sfi.GoToBegin()
	{
		if (image->GetPixel(it.GetIndex()))
		{
			variance += pow((image->GetPixel(it.GetIndex()) - mean), 2); //*(image->GetPixel(it.GetIndex()) - mean);
		}
	}
	sigma = sqrt(variance / numInteriorPixels1);
	
	float* vec = 0;
	vec = new float[2];
	vec[0] = mean;
	vec[1] = sigma;
	return vec;
}
static bool is_file_exist(const char *fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}
static OutputImageType::Pointer ThresholdFilterGT(InternalImageType::Pointer image)
{
	typedef itk::BinaryThresholdImageFilter<InternalImageType, OutputImageType> ThresholdType;
	ThresholdType::Pointer thresholdFilter = ThresholdType::New();
	//Watershed Segmentation Output as Threshold Output
	thresholdFilter->SetInput(image);
	thresholdFilter->SetInPlace(true);
	//Set Threshold filter limits
	//float value = image->GetPixel(seed);
	//std::cout << value << std::endl;
	thresholdFilter->SetLowerThreshold(1);
	//thresholdFilter->SetUpperThreshold(255);
	//Set Output Values
	thresholdFilter->SetOutsideValue(0);
	thresholdFilter->SetInsideValue(1);
	//  std::cout<<"Extracting Seed Label..."<<std::endl;
	t.Start();
	//Apply Threshold
	thresholdFilter->Update();
	t.Stop();
	//myfile2 << "...threshold Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	//myfile2 << "...threshold Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	// Compute some statistics
	return thresholdFilter->GetOutput();
}



#endif