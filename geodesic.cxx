#include "itkFastMarchingImageFilter.h"

#include "itkGeodesicActiveContourLevelSetImageFilter.h"

#include "utility.h"




static InternalImageType ::Pointer FastMarchingGeo(InternalImageType ::Pointer image, IndexType seed, char *filename)
{
	typedef  itk::FastMarchingImageFilter<InternalImageType, InternalImageType>    FastMarchingFilterType;
	FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
	
	fastMarching->SetSpeedConstant(1.0);

	fastMarching->SetOutputSize(image->GetBufferedRegion().GetSize());
	fastMarching->SetOutputSpacing(image->GetSpacing());
	fastMarching->SetOutputOrigin(image->GetOrigin());
	fastMarching->SetOutputDirection(image->GetDirection());

	typedef FastMarchingFilterType::NodeContainer NodeContainer;
	typedef FastMarchingFilterType::NodeType  NodeType;
	NodeContainer::Pointer seeds = NodeContainer::New();
	NodeType node;
	const double seedValue = -4.5;
	node.SetValue(seedValue);
	node.SetIndex(seed);
	seeds->Initialize();
	seeds->InsertElement(0, node);
	fastMarching->SetTrialPoints(seeds);
	
	myfile2 << "Performing FastMarching Segmentation..." << std::endl;
	//Perform Segmentation
	try
	{
		t.Start();
		fastMarching->Update();
		t.Stop();
	}
	catch (itk::ExceptionObject & excep)
	{
		std::cerr << "Exception caught in fast marching !" << std::endl;
		std::cerr << excep << std::endl;
	}
	//std::cout << "...FastMarching segmentation Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	myfile2 << "...FastMarching segmentation Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;

	myfile2 <<"GetSpacing " << fastMarching->GetOutput()->GetSpacing() << std::endl;
	myfile2<< "GetSpeedConstant "<< fastMarching->GetSpeedConstant() << std::endl;
	myfile2<<"GetStoppingValue "<< fastMarching->GetStoppingValue() << std::endl;

	string temp = "_FastMarchingGeo.nii.gz";
	string file = filename + temp;
	Write(fastMarching->GetOutput(), file);

	return fastMarching->GetOutput();
}
static InternalImageType::Pointer geodesicFilter(InternalImageType:: Pointer fimage, InternalImageType :: Pointer edgeImage, IndexType seed, 
														char *filename)
{
	typedef  itk::GeodesicActiveContourLevelSetImageFilter<InternalImageType, InternalImageType, InternalPixelType >   
																											GeodesicActiveContourFilterType;
		GeodesicActiveContourFilterType::Pointer geodesicActiveContour =GeodesicActiveContourFilterType::New();

		geodesicActiveContour->SetPropagationScaling(0.2);
		geodesicActiveContour->SetCurvatureScaling(1.0);
		geodesicActiveContour->SetAdvectionScaling(0.8);
		geodesicActiveContour->SetMaximumRMSError(0.02);
		geodesicActiveContour->SetNumberOfIterations(800);
		//geodesicActiveContour->SetReverseExpansionDirection(TRUE);
	//	geodesicActiveContour->SetAutoGenerateSpeedAdvection(TRUE);
		geodesicActiveContour->SetInput(fimage);		
		geodesicActiveContour->SetFeatureImage(edgeImage);
		
		myfile2 <<"staring geodesic"<< std::endl;
		try
		{
			t.Start();
			geodesicActiveContour->Update();
			t.Stop();
		}
		catch (itk::ExceptionObject & excep)
		{
			std::cerr << "Exception caught in geodesic !" << std::endl;
			std::cerr << excep << std::endl;
		}
		//std::cout << "...geodesic segmentation Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
		myfile2 << "...geodesic segmentation Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;

		std::cout << std::endl;

		myfile2 << "Max. no. iterations: " << geodesicActiveContour->GetNumberOfIterations() << std::endl;
		myfile2 << "Max. RMS error: " << geodesicActiveContour->GetMaximumRMSError() << std::endl;
		myfile2 << std::endl;
		myfile2 << "No. elpased iterations: " << geodesicActiveContour->GetElapsedIterations() << std::endl;
		myfile2 << "RMS change: " << geodesicActiveContour->GetRMSChange() << std::endl;
		myfile2 << "propogation: " << geodesicActiveContour->GetPropagationScaling() << std::endl;
		myfile2 << "Advection: " << geodesicActiveContour->GetAdvectionScaling() << std::endl;

		cout << "Max. no. iterations: " << geodesicActiveContour->GetNumberOfIterations() << std::endl;
		cout << "Max. RMS error: " << geodesicActiveContour->GetMaximumRMSError() << std::endl;
		cout << std::endl;
		cout << "No. elpased iterations: " << geodesicActiveContour->GetElapsedIterations() << std::endl;
		cout << "RMS change: " << geodesicActiveContour->GetRMSChange() << std::endl;
		//InternalImageType  speedImage;
		//speedImage = 
		//Write(geodesicActiveContour->GetSpeedImage(), );
		float value = geodesicActiveContour->GetOutput()->GetPixel(seed);
		std::cout << value << std::endl;	
		
		return geodesicActiveContour->GetOutput();
}

static OutputImageType::Pointer ThresholdFilterGeo(InternalImageType::Pointer image, IndexType seed)
{
	typedef itk::BinaryThresholdImageFilter<InternalImageType, OutputImageType> ThresholdType;
	ThresholdType::Pointer thresholdFilter = ThresholdType::New();
	//Watershed Segmentation Output as Threshold Output
	thresholdFilter->SetInput(image);
	thresholdFilter->SetInPlace(true);
	//Set Threshold filter limits
	float value = image->GetPixel(seed);
	std::cout << value << std::endl;
	//thresholdFilter->SetLowerThreshold(-1000.0);
	thresholdFilter->SetUpperThreshold(0);
	//Set Output Values
	thresholdFilter->SetOutsideValue(0);
	thresholdFilter->SetInsideValue(1);
	//  std::cout<<"Extracting Seed Label..."<<std::endl;
	t.Start();
	//Apply Threshold
	thresholdFilter->Update();
	t.Stop();
	//std::cout << "...threshold Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	myfile2 << "...threshold Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	// Compute some statistics
	return thresholdFilter->GetOutput();
}


