#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include <itkIndex.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include <itkNiftiImageIO.h>
#include "itkRescaleIntensityImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkWatershedImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkEllipsoidInteriorExteriorSpatialFunction.h"
#include "itkFloodFilledSpatialFunctionConditionalIterator.h"
#include "itkMaskImageFilter.h"
#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip> 
#include "../albertLib/toolsITK.h"
#include "../albertLib/toolsITK.cpp"
#include "itkTimeProbe.h"

//Template for the Dimension
//Image Properties
const unsigned int  Dimension = 3;
typedef float InternalPixelType;
typedef itk::Image<InternalPixelType, Dimension> InternalImageType;
typedef unsigned short OutputPixelType;
//typedef unsigned char OutputPixelType;
typedef itk::Image<OutputPixelType, Dimension> OutputImageType;
typedef  InternalImageType::IndexType  IndexType;
typedef  InternalImageType::IndexType  Indexarray[];
typedef long unsigned int WatershedPixelType;
typedef itk::Image<WatershedPixelType, Dimension> WatershedImageType;
class Segmentation
{
	itk::TimeProbe t;
	typedef itk::ImageFileReader<InternalImageType> ReaderType;  //the diacom image
	ReaderType::Pointer readImage;
public:
	InternalImageType::Pointer imageSeg;
	InternalImageType::RegionType region;
	InternalImageType::SizeType size;
	//constructor takes input and output path of image
public: Segmentation()
{         //instantiate reader
	readImage = ReaderType::New();
}
public: InternalImageType::Pointer ReadImage(char *filename)
{
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
	size = region.GetSize();
	std::cout << "Input image size: " << size << std::endl;
	return readImage->GetOutput();
}
public: int WriteImage(InternalImageType::Pointer image, std::string parameter)
{
	typedef itk::CastImageFilter< InternalImageType, OutputImageType > CastingFilterType1;
	CastingFilterType1::Pointer Castingfilter1 = CastingFilterType1::New();
	Castingfilter1->SetInput(image);
	cout << "Writting Image:: " << parameter << endl;
	try
	{
		toolsITK::writeDCMITKImage(Castingfilter1->GetOutput(), parameter);
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << e << std::endl;
		return 0;
	}
	return 1;
}
public: float getVolume(OutputImageType::Pointer imageStats, double sourceImageSpacing[])
{
	// Compute some statistics
	typedef itk::StatisticsImageFilter<OutputImageType> StatisticsImageFilterType;
	StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
	statisticsImageFilter->SetInput(imageStats);
	statisticsImageFilter->Update();
	float volume = (statisticsImageFilter->GetSum())*sourceImageSpacing[0] * sourceImageSpacing[1] * sourceImageSpacing[2];
	// Print statistics
	std::cout << std::endl;
	std::cout << "----------------------------------------------------" << std::endl;
	//std::cout << "Total number of voxels within lesion: " << statisticsImageFilter->GetSum() << std::endl;
	std::cout << "Total volume of segmented lesion: " << volume << " mm³" << std::endl;
	std::cout << "----------------------------------------------------" << std::endl;
	return volume;
}
public: float getVolume(InternalImageType::Pointer imageStats, double sourceImageSpacing[])
{
	// Compute some statistics
	typedef itk::StatisticsImageFilter<InternalImageType> StatisticsImageFilterType;
	StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
	statisticsImageFilter->SetInput(imageStats);
	statisticsImageFilter->Update();
	float volume = (statisticsImageFilter->GetSum())*sourceImageSpacing[0] * sourceImageSpacing[1] * sourceImageSpacing[2];
	// Print statistics
	std::cout << std::endl;
	std::cout << "----------------------------------------------------" << std::endl;
	//std::cout << "Total number of voxels within lesion: " << statisticsImageFilter->GetSum() << std::endl;
	std::cout << "Total volume of groundTruth: " << volume << " mm³" << std::endl;
	std::cout << "----------------------------------------------------" << std::endl;
	return volume;
}
public: InternalImageType::Pointer anisotropicDiffusion(InternalImageType::Pointer imagei, int iterations, float conduct)
{
	//Anisotropic Difussion (smooth the image)
	typedef itk::GradientAnisotropicDiffusionImageFilter<InternalImageType, InternalImageType> FilterType;
	FilterType::Pointer anisotropic = FilterType::New();
	//Input Image as Anisotropic Diffusion Input
	anisotropic->SetInput(imagei);
	anisotropic->SetInPlace(true);
	//Number of Iterations
	anisotropic->SetNumberOfIterations(iterations);
	//Step Time
	//must be smaller than  = imagespacing / 2^(dimention +1) used is step  = imagespacing / 2^(dimention +2)
	float step = imagei->GetSpacing()[0] / pow(2, Dimension + 2);
	//std::cout<< "the step is :: " << step ;
	anisotropic->SetTimeStep(step);
	//Conductance
	anisotropic->SetConductanceParameter(conduct);
	std::cout << "Applying Anistropic Diffusion ..." << std::endl;
	//Apply Anisotropic Difussion
	t.Start();
	anisotropic->Update();
	t.Stop();
	std::cout << "...Anisotropic diffusion Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	return anisotropic->GetOutput();
}

public: InternalImageType::Pointer normalizeGtImage(InternalImageType::Pointer image)
{
	typedef itk::RescaleIntensityImageFilter<InternalImageType, InternalImageType> RescaleFilterType;
	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	//Image GrayValue Limits
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(1);
	//Anisotropic Diffusion Output as Rescaler Input
	rescaler->SetInput(image);
	rescaler->SetInPlace(true);
	//  std::cout<<"Rescaling Image..."<<std::endl;
	//Rescale the Image
	t.Start();
	rescaler->Update();
	t.Stop();
	std::cout << "...Rescale Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	return rescaler->GetOutput();
}
public: InternalImageType::Pointer Smoothing(InternalImageType::Pointer image, float sigma)
{   //Gradient Magnitude Filter (generates the height function)
	typedef itk::RecursiveGaussianImageFilter<InternalImageType, InternalImageType>smoothingFilterType;
	smoothingFilterType::Pointer smoothingFilterX = smoothingFilterType::New();
	smoothingFilterType::Pointer smoothingFilterY = smoothingFilterType::New();
	smoothingFilterType::Pointer smoothingFilterZ = smoothingFilterType::New();

	smoothingFilterX->SetDirection(0);    //0 --> X direction
	smoothingFilterY->SetDirection(1);    //1 --> Y direction
	smoothingFilterZ->SetDirection(2);    //2 --> Z direction
	smoothingFilterX->SetOrder(smoothingFilterType::ZeroOrder);
	smoothingFilterY->SetOrder(smoothingFilterType::ZeroOrder);
	smoothingFilterZ->SetOrder(smoothingFilterType::ZeroOrder);

	//Rescaler Output as  Filter Input
	smoothingFilterX->SetNormalizeAcrossScale(true);
	smoothingFilterY->SetNormalizeAcrossScale(true);
	smoothingFilterZ->SetNormalizeAcrossScale(true);

	smoothingFilterX->SetInput(image);
	smoothingFilterY->SetInput(smoothingFilterX->GetOutput());
	smoothingFilterZ->SetInput(smoothingFilterY->GetOutput());


	//Sigma (standard deviation of Gaussian)
	smoothingFilterX->SetSigma(sigma);  ///parameter to optimize
	smoothingFilterY->SetSigma(sigma);  ///parameter to optimize
	smoothingFilterZ->SetSigma(sigma);  ///parameter to optimize

	std::cout << "Applying gaussian smoothing Filter..." << sigma << std::endl;
	//Apply Gradient Magnitude Filter
	try
	{
		t.Start();
		smoothingFilterZ->Update();
		t.Stop();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e << std::endl;
		return NULL;
	}
	std::cout << "...Smoothing Filter Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	return smoothingFilterZ->GetOutput();
}
public: InternalImageType::Pointer HeightFunction(InternalImageType::Pointer image)
{
	//Gradient Magnitude Filter (generates the height function)
	typedef itk::GradientMagnitudeImageFilter<InternalImageType, InternalImageType>GradientFilterType;
	GradientFilterType::Pointer gradienFilter = GradientFilterType::New();
	gradienFilter->SetInput(image);
	//  gradienFilter->SetSigma( sigma ); 
	// gradienFilter->SetNormalizeAcrossScale(false);
	//Apply Gradient Magnitude Filter
	std::cout << "Applying Gradient Magnitude Filter ..." << std::endl;
	try
	{
		t.Start();
		gradienFilter->Update();
		t.Stop();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e << std::endl;
		return NULL;
	}
	std::cout << "...Gradient Magnitude Filter Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	return gradienFilter->GetOutput();
}
public: InternalImageType::Pointer CreateMask(InternalImageType::Pointer image, InternalImageType::Pointer maskingImage, IndexType axis, IndexType seedPos)
{
	// Image size and spacing parameters
	unsigned long xExtent = size[0];
	unsigned long yExtent = size[1];
	unsigned long zExtent = size[2];
	unsigned long sourceImageSize[] = { xExtent, yExtent, zExtent };
	double sourceImageSpacing[] = { image->GetSpacing()[0],image->GetSpacing()[1],image->GetSpacing()[2] };  // (pixels)
																											 //std::cout<< "sourceImageSpacing: " << sourceImageSpacing[0] << " x " << sourceImageSpacing[1] << " x " << sourceImageSpacing[2] << " mm^3" << std::endl;
	double sourceImageOrigin[] = { image->GetOrigin()[0],image->GetOrigin()[1],image->GetOrigin()[2] };  // (mm/pixels)
																										 //std::cout<< "sourceImageOrigin: " << sourceImageOrigin[0] << " x " << sourceImageOrigin[1] << " x " << sourceImageOrigin[2] << " mm^3" << std::endl;
																										 // Creates the sourceImage (but doesn't set the size or allocate memory)
	InternalImageType::Pointer sourceImage = InternalImageType::New();
	sourceImage->SetOrigin(sourceImageOrigin);
	sourceImage->SetSpacing(sourceImageSpacing);
	//std::cout<< "sourceImageSpacing: " << sourceImageSpacing[0] << " x " << sourceImageSpacing[1] << " x " << sourceImageSpacing[2] << " mm^3" << std::endl;
	//-----The following block allocates the sourceImage-----

	InternalImageType::SizeType sourceImageSizeObject;   // Create a size object native to the sourceImage type
	InternalImageType::RegionType largestPossibleRegion; // Create a region object native to the sourceImage type

	sourceImageSizeObject.SetSize(sourceImageSize);             // Set the size object to the array defined earlier
	largestPossibleRegion.SetSize(sourceImageSizeObject);       // Resize the region
	sourceImage->SetLargestPossibleRegion(largestPossibleRegion); // Set the largest legal region size (i.e. the size of the whole sourceImage) to what we just defined
	sourceImage->SetBufferedRegion(largestPossibleRegion);      // Set the buffered region
	sourceImage->SetRequestedRegion(largestPossibleRegion);     // Set the requested region
	sourceImage->Allocate();                                      // Now allocate memory for the sourceImage

																  //std::cout << "New physical sourceImage allocated\n";
																  //################################ Create ellipsoid in sourceImage ##################################
																  // Ellipsoid spatial function typedef
	typedef itk::EllipsoidInteriorExteriorSpatialFunction<Dimension> EllipsoidFunctionType;
	EllipsoidFunctionType::Pointer spatialFunc = EllipsoidFunctionType::New();// Create an ellipsoid spatial function for the source image

																			  // Define and set the axes lengths for the ellipsoid
	EllipsoidFunctionType::InputType axes;
	axes[0] = axis[0];             //mm
	axes[1] = axis[1];             //mm
	axes[2] = axis[2];             //mm
	spatialFunc->SetAxes(axes);

	//ReaderWriter/ Define and set the center of the ellipsoid in physical space
	EllipsoidFunctionType::InputType center;
	center[0] = seedPos[0] * sourceImageSpacing[0];  // world coordinates(mm) -> pixels * spacing
	center[1] = seedPos[1] * sourceImageSpacing[1];  // world coordinates(mm) -> pixels * spacing
	center[2] = seedPos[2] * sourceImageSpacing[2];  // world coordinates(mm) -> pixels * spacing
	spatialFunc->SetCenter(center);

	double data[] = { 0, 1, 0, 1, 0, 0, 0, 0, 1 };
	vnl_matrix<double> orientations(data, 3, 3);

	// Set the orientations of the ellipsoids
	spatialFunc->SetOrientations(orientations);
	typedef  IndexType::IndexValueType   IndexValueType;
	IndexType seed;
	seed[0] = static_cast<IndexValueType>(seedPos[0]); // (pixels)
	seed[1] = static_cast<IndexValueType>(seedPos[1]); // (pixels)
	seed[2] = static_cast<IndexValueType>(seedPos[2]); // (pixels)

	itk::FloodFilledSpatialFunctionConditionalIterator<InternalImageType, EllipsoidFunctionType>
		sfi = itk::FloodFilledSpatialFunctionConditionalIterator<InternalImageType,
		EllipsoidFunctionType>(sourceImage, spatialFunc, seed);

	std::cout << "Seed position : " << seed << std::endl;

	// Iterate through the entire image and set interior pixels to 255
	int numInteriorPixels1 = 0;
	unsigned char interiorPixelValue = 1;
	// double intensity = 0;
	for (sfi.GoToBegin(); !sfi.IsAtEnd(); ++sfi) // sfi.GoToBegin()
	{
		//std::cout << "Index evaluated: " << sfi.GetIndex() << std::endl;
		sfi.Set(interiorPixelValue);
		++numInteriorPixels1;
		//intensity += maskingImage->GetPixel(sfi.GetIndex());
	}
	cout << "Returning the created image" << endl;
	return getGTMask(maskingImage, sourceImage);
}

public: OutputImageType::Pointer watershedSegmentation(InternalImageType::Pointer image, float threshold, float level, IndexType seed, char *filename, int iterations)
{
	//Watershed segmentation
	typedef itk::WatershedImageFilter<InternalImageType> WatershedFilterType;
	WatershedFilterType::Pointer watershedFilter = WatershedFilterType::New();
	//Gradient Magnitude Filter Output as Watershed Input
	watershedFilter->SetInput(image);
	//Threshold
	watershedFilter->SetThreshold(threshold);
	//Level
	watershedFilter->SetLevel(level);
	std::cout << "Performing Watershed Segmentation..." << std::endl;
	//Perform Segmentation
	try
	{
		t.Start();
		watershedFilter->Update();
		t.Stop();
	}
	catch (itk::ExceptionObject & excep)
	{
		std::cerr << "Exception caught here !" << std::endl;
		std::cerr << excep << std::endl;
	}
	OutputPixelType value = watershedFilter->GetOutput()->GetPixel(seed);
	std::cout << "the tresholded intensity is " << value << endl;
	float newlevel = watershedFilter->GetLevel();;
	if (value == 1)
	{
		while (value <= 1 && newlevel > 0)
		{
			newlevel = newlevel - 0.1;
			std::cout << "new level value" << newlevel << endl;
			std::cout << "Performing Watershed Segmentation..." << std::endl;
			watershedFilter->SetLevel(newlevel);
			try {
				t.Start();
				watershedFilter->Update();
				t.Stop();
			}
			catch (itk::ExceptionObject & excep)
			{
				std::cerr << "Exception caught here !" << std::endl;
				std::cerr << excep << std::endl;
			}

			value = watershedFilter->GetOutput()->GetPixel(seed);
			std::cout << "the tresholded intensity is " << value << endl;
			newlevel = watershedFilter->GetLevel();
		}
		newlevel = newlevel + 0.1;
		value = 1;
		while (value <= 1 && newlevel > 0)
		{
			newlevel = newlevel - 0.05;
			std::cout << "new level value" << newlevel << endl;
			std::cout << "Performing Watershed Segmentation..." << std::endl;
			watershedFilter->SetLevel(newlevel);
			try {
				t.Start();
				watershedFilter->Update();
				t.Stop();
			}
			catch (itk::ExceptionObject & excep)
			{
				std::cerr << "Exception caught here !" << std::endl;
				std::cerr << excep << std::endl;
			}

			value = watershedFilter->GetOutput()->GetPixel(seed);
			std::cout << "the tresholded intensity is " << value << endl;
			newlevel = watershedFilter->GetLevel();
		}

		newlevel = newlevel + 0.05;
		value = 1;
		while (value <= 1 && newlevel > 0)
		{
			newlevel = newlevel - 0.01;
			std::cout << "new level value" << newlevel << endl;
			std::cout << "Performing Watershed Segmentation..." << std::endl;
			watershedFilter->SetLevel(newlevel);
			try {
				t.Start();
				watershedFilter->Update();
				t.Stop();
			}
			catch (itk::ExceptionObject & excep)
			{
				std::cerr << "Exception caught here !" << std::endl;
				std::cerr << excep << std::endl;
			}

			value = watershedFilter->GetOutput()->GetPixel(seed);
			std::cout << "the tresholded intensity is " << value << endl;
			newlevel = watershedFilter->GetLevel();
		}
		std::cout << "the tresholded intensity is " << value << endl;
	}
	std::cout << "...Watershed segmentation Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	std::stringstream ss;
	ss << "_"<<iterations<<"_"<< newlevel<<"_"<< value;
	std::string temp = ss.str();
	string waterShedfile = filename + temp + "watershed.dcm";
	//WriteWatershedIntermediateimage(watershedFilter->GetOutput(),waterShedfile);
	typedef itk::CastImageFilter< WatershedImageType, OutputImageType > CastingFilterType2;
	CastingFilterType2::Pointer Castingfilter2 = CastingFilterType2::New();
	Castingfilter2->SetInput(watershedFilter->GetOutput());
	Castingfilter2->Update();
	//Seed Label
	// WriteDiacomImage(watershedFilter->GetOutput(),waterShedfile);
	try
	{
		toolsITK::writeDCMITKImage(Castingfilter2->GetOutput(), waterShedfile);
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file writer " << std::endl;
		std::cerr << e << std::endl;
		return 0;
	}
	OutputPixelType valueInten = Castingfilter2->GetOutput()->GetPixel(seed);
	cout << "my value is   " << valueInten << endl;
	//Label Threshold
	typedef itk::BinaryThresholdImageFilter<OutputImageType, OutputImageType> ThresholdType;
	ThresholdType::Pointer thresholdFilter = ThresholdType::New();
	//Watershed Segmentation Output as Threshold Output
	thresholdFilter->SetInput(Castingfilter2->GetOutput());
	thresholdFilter->SetInPlace(true);
	//Set Threshold filter limits
	thresholdFilter->SetLowerThreshold(valueInten);
	thresholdFilter->SetUpperThreshold(valueInten);
	//Set Output Values
	thresholdFilter->SetOutsideValue(0);
	thresholdFilter->SetInsideValue(1);
	//  std::cout<<"Extracting Seed Label..."<<std::endl;
	t.Start();
	//Apply Threshold
	thresholdFilter->Update();
	t.Stop();
	// Compute some statistics
	double sourceImageSpacing[] = { image->GetSpacing()[0],image->GetSpacing()[1],image->GetSpacing()[2] };  // (pixels)
																											 //std::cout<< "Image spacing: " << sourceImageSpacing[0] << " x " << sourceImageSpacing[1] << " x " << sourceImageSpacing[2] << " mm^3" << std::endl;
	getVolume(thresholdFilter->GetOutput(), sourceImageSpacing);
	ss.str(""); temp.clear(); waterShedfile.clear();
	ss << "_" << iterations; temp = ss.str();
	waterShedfile = filename + temp + "_segmented.dcm";
	try
	{
		toolsITK::writeDCMITKImage(thresholdFilter->GetOutput(), waterShedfile);
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file writer " << std::endl;
		std::cerr << e << std::endl;
		return NULL;
	}
	return thresholdFilter->GetOutput();
}

public: int** getCoordinates(IndexType seed, IndexType axex, double sourceImageSpacing[])
{
	int** p = 0; p = new int*[8];
	for (int h = 0; h < 8; h++)
	{
		p[h] = new int[3];
	}

	axex[0] = axex[0] / (2 * sourceImageSpacing[0]);          //pixels as axes is in mm
	axex[1] = axex[1] / (2 * sourceImageSpacing[1]);
	axex[2] = axex[2] / (2 * sourceImageSpacing[2]);

	//cout<< axex <<endl;
	// cout<< seed <<endl;
	p[0][0] = seed[0] - axex[0];                          //pixels as axes is in mm
	p[0][1] = seed[1] - axex[1];
	p[0][2] = seed[2] - axex[2];

	p[1][0] = seed[0] + axex[0];
	p[1][1] = seed[1] - axex[1];
	p[1][2] = seed[2] - axex[2];

	p[2][0] = seed[0] + axex[0];
	p[2][1] = seed[1] + axex[1];
	p[2][2] = seed[2] - axex[2];

	p[3][0] = seed[0] - axex[0];
	p[3][1] = seed[1] + axex[1];
	p[3][2] = seed[2] - axex[2];

	p[4][0] = seed[0] - axex[0];
	p[4][1] = seed[1] - axex[1];
	p[4][2] = seed[2] + axex[2];

	p[5][0] = seed[0] + axex[0];
	p[5][1] = seed[1] - axex[1];
	p[5][2] = seed[2] + axex[2];

	p[6][0] = seed[0] + axex[0];
	p[6][1] = seed[1] + axex[1];
	p[6][2] = seed[2] + axex[2];

	p[7][0] = seed[0] - axex[0];
	p[7][1] = seed[1] + axex[1];
	p[7][2] = seed[2] + axex[2];

	for (int i = 0; i<8; i++)
	{
		for (int j = 0; j <3; j++)
		{
			if (p[i][j] <= 1)
				p[i][j] = axex[j];
			if (p[i][j] >= size[j])
				p[i][j] = size[j] - axex[j];
		}
	}
	return p;
}

public: InternalImageType::Pointer getGTMask(InternalImageType::Pointer image, InternalImageType::Pointer maskingImage)
{
	typedef itk::MaskImageFilter< InternalImageType, InternalImageType > MaskFilterType;
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();

	maskFilter->SetInput(image);
	maskFilter->SetMaskImage(maskingImage);

	maskFilter->Update();
	return maskFilter->GetOutput();
}

public: float* getCnrSNR(InternalImageType::Pointer image)
{
	float intensity = 0.0;
	int numInteriorPixels1 = 0;
	//InternalImageType::Pointer im = maskingImage;
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
			variance += (image->GetPixel(it.GetIndex()) - mean)*(image->GetPixel(it.GetIndex()) - mean);
		}
	}
	// std::cout << "variance of the region is "<<variance << std::endl;
	sigma = sqrt(variance / numInteriorPixels1);
	// std::cout << "Sigma of the region is "<<sigma << std::endl;
	float* vec = 0;
	vec = new float[2];
	vec[0] = mean;
	vec[1] = sigma;
	return vec;
}
		bool is_file_exist(const char *fileName)
		{
			std::ifstream infile(fileName);
			return infile.good();
		}
public: InternalImageType::Pointer MedianSmoothing(InternalImageType::Pointer iImage, int iRadius)
{   //initialize Median Filter
	typedef itk::MedianImageFilter<InternalImageType, InternalImageType > MedianFilterType;
	MedianFilterType::Pointer medianFilter = MedianFilterType::New();
	//create radius
	MedianFilterType::InputSizeType radius;
	radius.Fill(iRadius);
	//set radius within filter
	medianFilter->SetRadius(radius);
	//set input Image
	medianFilter->SetInput(iImage);
	try
	{
		medianFilter->Update();
	}
	catch (itk::ExceptionObject & e) {
		std::cerr << "Median Filter Error" << std::endl;
		std::cerr << e << std::endl;
	}
	return medianFilter->GetOutput();
}
public: float* evaluate(InternalImageType::Pointer imageGt, OutputImageType::Pointer imageSeg)
{
	//cout<<"11"<<endl;
	typedef itk::CastImageFilter< InternalImageType, WatershedImageType > CastingFilterType1;
	CastingFilterType1::Pointer Castingfilter1 = CastingFilterType1::New();
	Castingfilter1->SetInput(imageGt);
	Castingfilter1->Update();

	typedef itk::CastImageFilter< OutputImageType, WatershedImageType > CastingFilterType2;
	CastingFilterType2::Pointer Castingfilter2 = CastingFilterType2::New();
	Castingfilter2->SetInput(imageSeg);
	Castingfilter2->Update();

	typedef itk::LabelOverlapMeasuresImageFilter<WatershedImageType> FilterType;
	FilterType::Pointer labelFilter = FilterType::New();
	labelFilter->SetSourceImage(Castingfilter1->GetOutput());
	labelFilter->SetTargetImage(Castingfilter2->GetOutput());
	try
	{
		labelFilter->Update();
	}
	catch (itk::ExceptionObject & e) {
		std::cerr << "label Filter Error" << std::endl;
		std::cerr << e << std::endl;
	}
	std::cout << std::setw(17) << "Union (jaccard)" << std::setw(17) << labelFilter->GetJaccardCoefficient() << std::endl;
	std::cout << std::setw(17) << "A^B" << std::setw(17) << labelFilter->GetTotalOverlap() << std::endl;
	std::cout << std::setw(17) << "Mean (dice)" << std::setw(17) << labelFilter->GetDiceCoefficient() << std::endl;
	std::cout << std::setw(17) << "Volume sim." << std::setw(17) << labelFilter->GetVolumeSimilarity() << std::endl;
	std::cout << std::setw(17) << "False negative" << std::setw(17) << labelFilter->GetFalseNegativeError() << std::endl;
	std::cout << std::setw(17) << "False positive" << std::setw(17) << labelFilter->GetFalsePositiveError() << std::endl;
	float* vec = 0;
	vec = new float[5];
	vec[0] = labelFilter->GetJaccardCoefficient();
	vec[1] = labelFilter->GetDiceCoefficient();
	vec[2] = labelFilter->GetFalseNegativeError();
	vec[3] = labelFilter->GetFalsePositiveError();
	vec[4] = labelFilter->GetVolumeSimilarity();
	return vec;
}

};

int main(int argc, char *argv[])
{
	//cout<< argv[10]<<endl;
	if (argc != 9) {
		std::cout << "Input_image Output_image SeedX SeedY SeedZ" << std::endl;
		return 0;
	}


	itk::TimeProbe clock;
	clock.Start();
	Segmentation watershedseg;
	Segmentation tempObj;
	char * name = argv[6];
	char * file = argv[8];
	std::string newName;

	int nei = 5;                    //mm diameter for each background  pixel     
	float conductance = 3;
	float level = 0.80;               // Controls watershed depth. Percentage(0.0-1.0) of the maximum depth in the input image.
	int iterations = 2;              // Typical 5. more interactions will smooth further but will increase computing time.
	float sigma = 0.1;
	float threshold = 0.005;        // 0.1 - controls the lower thresholding of the input. Percentage(0.0-1.0) of the maximum depth in the input image.
									// (A rule of thumb is to set the Threshold to be about 1 / 100 of the Level.)
	IndexType seedPos;
	IndexType seedPoint;
	seedPoint[0] = atoi(argv[2]); // (pixels)
	seedPoint[1] = atoi(argv[3]); // (pixels)
	seedPoint[2] = atoi(argv[4]); // (pixels)

	IndexType axesImage;
	axesImage[0] = atoi(argv[5]);      // (mm)
	axesImage[1] = atoi(argv[5]);      // (mm)
	axesImage[2] = atoi(argv[5]);      // (mm)

	IndexType axes;
	axes[0] = 30;      // (mm)
	axes[1] = 30;      // (mm)
	axes[2] = 30;      // (mm)

					   //background mean and stnadard deviation
	float ** bgMuStd = 0;
	bgMuStd = new float*[8];
	for (int h = 0; h < 8; h++) { bgMuStd[h] = new float[2]; }
	float bgMean = 0, bgSigma = 0;
	//ground truth mean and standard deviation
	float* sigMuStd = 0;
	sigMuStd = new float[2];

	//to store dice and other evaluation measure
	float *labelMeasure = 0;
	labelMeasure = new float[5];

	IndexType neighbours;
	neighbours[0] = nei;            // (pixels)
	neighbours[1] = nei;            // (pixels)
	neighbours[2] = nei;            // (pixels)

									//Perform algorithm
	InternalImageType::Pointer  backGroundImages, diacomImage, gtImage, gtImageN, maskImageNew, medImgFilter, medImgGtFilter, smoothimage, gtSmoothImage;
	OutputImageType::Pointer segmentationImage;
	std::stringstream sss;
	std::string stn;

	//load the image 

	diacomImage = watershedseg.ReadImage(argv[1]);
	medImgFilter = watershedseg.MedianSmoothing(diacomImage, 1);  //added to remove the salt and pepper noise radius 1 means 3*3  mask
																  //get Image Spacing
	double sourceImageSpacing[] = { diacomImage->GetSpacing()[0],diacomImage->GetSpacing()[1],diacomImage->GetSpacing()[2] }; //get the spacing of the image

																															  //load the ground truth
	gtImage = tempObj.ReadImage(argv[7]);
	gtImageN = tempObj.normalizeGtImage(gtImage);
	//creating the groundtruth Mask iamge to get the mean and standard deviation of the ground truth
	InternalImageType::Pointer gtMaskImage = watershedseg.getGTMask(medImgFilter, gtImageN);
	sss << "_"; stn = sss.str();        string hei = name + stn + "maskGT.dcm";
	tempObj.WriteImage(gtMaskImage, hei);

	//medImgGtFilter = tempObj.MedianSmoothing(gtImage, 2); 
	//Truncating the whole image to get a mask od 50mm
	sss << "Mask.dcm"; stn = sss.str(); newName = name + stn;
	char *tempname = new char[newName.length() + 1];
	std::strcpy(tempname, newName.c_str());
	cout << tempname << endl;
	maskImageNew = watershedseg.CreateMask(diacomImage, medImgFilter, axesImage, seedPoint);
	watershedseg.WriteImage(maskImageNew, newName);
	newName.clear(); sss.str(""); stn.clear();
	//create the background mask to get the mean and standard deviation of the background pixels to calculate CNR
	//generate the 8 co-ordinates around the mask for background information
	int** cubeCordinates = watershedseg.getCoordinates(seedPoint, axes, sourceImageSpacing);
	//std::cout<< "cubeCordinates:: "<< cubeCordinates[0][0]<<"  " <<cubeCordinates[0][1]<<"  " <<cubeCordinates[0][2]<<endl;
	for (int i = 0; i<8; i++)
	{
		seedPos[0] = cubeCordinates[i][0];
		seedPos[1] = cubeCordinates[i][1];
		seedPos[2] = cubeCordinates[i][2];
		sss << "_" << seedPos[0] << "_" << seedPos[1] << "_" << seedPos[2] << "_Mask.dcm";  stn = sss.str();   newName = name + stn;
		char *tempnameN = new char[newName.length() + 1];
		backGroundImages = watershedseg.CreateMask(diacomImage, maskImageNew, neighbours, seedPos);
		watershedseg.WriteImage(backGroundImages, newName);
		newName.clear(); sss.str(""); stn.clear();
		//  cout<<newName<<"\t"<<sss<<"\t"<<stn;

		bgMuStd[i] = watershedseg.getCnrSNR(backGroundImages);
	}

	for (int i = 0; i<8; i++)
	{
		bgMean += bgMuStd[i][0];
		bgSigma += bgMuStd[i][1];
	}
	bgMean = bgMean / 8;
	bgSigma = bgSigma / 8;


	sigMuStd = watershedseg.getCnrSNR(gtMaskImage);

	float cnr = abs(sigMuStd[0] - bgMean) / sqrt(sigMuStd[1] * sigMuStd[1] + bgSigma*bgSigma);  //ref in draft
	float snr = sigMuStd[0] / bgSigma;

	std::cout << "###############################------------" << sigMuStd[0] << " , " << sigMuStd[1] << ", " << bgMean << " , " << bgSigma << std::endl;
	std::cout << "###############################----CNR-----" << cnr << endl;
	std::cout << "###############################----SNR-----" << snr << endl;
	ofstream myfile2;
	myfile2.open(file, ios::out | ios::app);
	myfile2 << name << "," << "," << sigMuStd[0] << "," << sigMuStd[1] << "," << bgMean << " , " << bgSigma << "," << cnr << "," << snr << "\n";
	myfile2.close();
	for (int ii = 0; ii <20; ii++)
	{
		//Anisotropic diffusion on the image and gt
		sss << "_" << iterations << "_AD.dcm";     stn = sss.str();         newName = name + stn;
		char *tempname1 = new char[newName.length() + 1];
		std::strcpy(tempname1, newName.c_str());
		//cout<<tempname1<<endl;
		smoothimage = watershedseg.anisotropicDiffusion(maskImageNew, iterations, conductance);
		watershedseg.WriteImage(smoothimage, newName);

		newName.clear(); sss.str(""); stn.clear();

		// performing smoothing on the image and gt
		/*InternalImageType::Pointer gaussianSmoothimage, gaussianSmoothimageGt;
		sss << "_" << iterations << "_" << sigma << "_GS.dcm";       stn = sss.str();        newName = name + stn;
		char *tempname2 = new char[newName.length() + 1];
		std::strcpy(tempname2, newName.c_str());

		gaussianSmoothimage = watershedseg.Smoothing(smoothimage, sigma);
		watershedseg.WriteImage(gaussianSmoothimage, newName);

		newName.clear(); sss.str(""); stn.clear();*/

	

		gtMaskImage = watershedseg.getGTMask(smoothimage, gtImageN);
		sss << "_"<<iterations<<"_"; stn = sss.str();        string hei = name + stn + "maskGTADGS.dcm";
		tempObj.WriteImage(gtMaskImage, hei);

		for (int i = 0; i<8; i++)
		{
			seedPos[0] = cubeCordinates[i][0];
			seedPos[1] = cubeCordinates[i][1];
			seedPos[2] = cubeCordinates[i][2];
			sss << "_" << iterations << "_" << seedPos[0] << "_" << seedPos[1] << "_" << seedPos[2] << "_Mask.dcm";  stn = sss.str();   newName = name + stn;
			char *tempnameN = new char[newName.length() + 1];
			backGroundImages = watershedseg.CreateMask(diacomImage, smoothimage, neighbours, seedPos);
			watershedseg.WriteImage(backGroundImages, newName);
			newName.clear(); sss.str(""); stn.clear();
			//  cout<<newName<<"\t"<<sss<<"\t"<<stn;
			bgMuStd[i] = watershedseg.getCnrSNR(backGroundImages);
		}
		for (int i = 0; i<8; i++)
		{
			bgMean += bgMuStd[i][0];
			bgSigma += bgMuStd[i][1];
		}
		bgMean = bgMean / 8;
		bgSigma = bgSigma / 8;

		sigMuStd = watershedseg.getCnrSNR(gtMaskImage);

		float cnrnew = abs(sigMuStd[0] - bgMean) / sqrt(sigMuStd[1] * sigMuStd[1] + bgSigma*bgSigma);  //ref in draft
		float snrnew = sigMuStd[0] / bgSigma;
		std::cout << "###############################------------" << sigMuStd[0] << " , " << sigMuStd[1] << ", " << bgMean << " , " << bgSigma << std::endl;
		std::cout << "###############################----NewCnr-----" << cnrnew << endl;
		std::cout << "###############################----SNRSnr-----" << snrnew << endl;
		myfile2.open(file , ios::out | ios::app);
		myfile2 << name <<"," << iterations <<"," << sigMuStd[0] << "," << sigMuStd[1] << "," << bgMean << " , " << bgSigma << "," << cnrnew << "," << snrnew << ",";
		myfile2.close();
		
		// calculating the image and gt f¡gradient image or height function
		InternalImageType::Pointer  heightImage;
		sss << "_" << iterations << "_GM.dcm";     stn = sss.str();        newName = name + stn;
		char *temp = new char[newName.length() + 1];
		std::strcpy(temp, newName.c_str());

		heightImage = watershedseg.HeightFunction(smoothimage);
		watershedseg.WriteImage(heightImage, newName);

		newName.clear(); sss.str(""); stn.clear();

		// creating the mask again for the image to avoide the edges artifact due to masking
		InternalImageType::Pointer maskSmoothImage;
		sss << "_" << iterations << "_GM_Mask.dcm";        stn = sss.str();        newName = name + stn;
		char *tempname3 = new char[newName.length() + 1];
		std::strcpy(tempname3, newName.c_str());
		//cout<<tempname3<<endl;
		maskSmoothImage = watershedseg.CreateMask(diacomImage, heightImage, axes, seedPoint);
		watershedseg.WriteImage(maskSmoothImage, newName);

		newName.clear(); sss.str(""); stn.clear();
		segmentationImage = watershedseg.watershedSegmentation(maskSmoothImage, threshold, level, seedPoint,name,iterations);
		float volumeSeg = watershedseg.getVolume(segmentationImage, sourceImageSpacing);
		float volumeGt = watershedseg.getVolume(gtImage, sourceImageSpacing);
		labelMeasure = watershedseg.evaluate(gtImage, segmentationImage);
		myfile2.open(file , ios::out | ios::app);
		myfile2 << labelMeasure[0] << "," << labelMeasure[1] << "," << labelMeasure[2] << "," << labelMeasure[3] << "," << volumeGt << "," << volumeSeg << "," << labelMeasure[4] << "\n";
		myfile2.close();
		
		iterations = iterations + 2;
		maskImageNew = smoothimage;
		
		// gtImageN = gtSmoothImage; 

	}

	// print computation time at screen
	clock.Stop();
	std::cout << "###############################" << std::endl;
	std::cout << "SegmentationABUS" << std::endl;
	std::cout << "Total time: " << clock.GetTotal() / 60.0 << " min." << std::endl;
	std::cout << "Total time: " << clock.GetTotal() << " s." << std::endl;
	std::cout << "###############################" << std::endl;

	return 0;
}



