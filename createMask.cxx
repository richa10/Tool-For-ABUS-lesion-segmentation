

#include "itkEllipsoidInteriorExteriorSpatialFunction.h"
#include "itkFloodFilledSpatialFunctionConditionalIterator.h"
#include "itkMaskImageFilter.h"

#include "utility.h"

static InternalImageType::Pointer getMask(InternalImageType::Pointer image, InternalImageType::Pointer maskingImage)
{
	typedef itk::MaskImageFilter< InternalImageType, InternalImageType > MaskFilterType;
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();
	//cout << "Here";
	
	maskFilter->SetInput(image);
	//cout << "Here1";
	maskFilter->SetMaskImage(maskingImage);
	//cout << "Here3";
	try {
		maskFilter->Update();
	}
	catch (itk::ExceptionObject & excep)
	{
		std::cerr << "Exception caught in Get mask !" << std::endl;
		std::cerr << excep << std::endl;
	}
	return maskFilter->GetOutput();
}
static InternalImageType::Pointer getMask(InternalImageType::Pointer image, OutputImageType::Pointer maskingImage)
{
	typedef itk::MaskImageFilter< InternalImageType, OutputImageType > MaskFilterType;
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();

	maskFilter->SetInput(image);
	maskFilter->SetMaskImage(maskingImage);

	maskFilter->Update();
	return maskFilter->GetOutput();
}

static InternalImageType :: Pointer CreateMask(InternalImageType :: Pointer image, InternalImageType::Pointer maskImage ,IndexType axis, IndexType seedPos)
{
	// Image size and spacing parameters
	unsigned long xExtent = Imagesize[0];
	unsigned long yExtent = Imagesize[1];
	unsigned long zExtent = Imagesize[2];
	unsigned long sourceImageSize[] = { xExtent, yExtent, zExtent };
	double sourceImageSpacing[] = { image->GetSpacing()[0],image->GetSpacing()[1],image->GetSpacing()[2] };  // (pixels)
	//std::cout << "sourceImageSpacing: " << sourceImageSpacing[0] << " x " << sourceImageSpacing[1] << " x " << sourceImageSpacing[2] << " mm^3" << std::endl;

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

	sourceImage->Update();
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

	//std::cout << "Seed position : " << seed << std::endl;

	// Iterate through the entire image and set interior pixels to 1
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


	//cout << "Returning the mask image" << endl;

	return getMask(maskImage, sourceImage);
}

