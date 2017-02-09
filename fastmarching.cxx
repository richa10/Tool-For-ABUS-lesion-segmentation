#include "itkFastMarchingImageFilter.h"
#include "itkInvertIntensityImageFilter.h."
#include <itkMinimumMaximumImageCalculator.h>
#include "utility.h"

static InternalImageType::Pointer FastMarching(InternalImageType::Pointer image, IndexType seed, 
											char *filename, const double stoppingTime)
{
	
	typedef  itk::FastMarchingImageFilter<InternalImageType, InternalImageType>    FastMarchingFilterType;
	FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
	
	fastMarching->SetInput(image);
		
	typedef FastMarchingFilterType::NodeContainer NodeContainer;
	typedef FastMarchingFilterType::NodeType  NodeType;
	NodeContainer::Pointer seeds = NodeContainer::New();
	NodeType node;
	const double seedValue = 0.0;

	node.SetValue(seedValue);
	node.SetIndex(seed);

	seeds->Initialize();
	seeds->InsertElement(0, node);
	//seeds->InsertElement(1, node);

	fastMarching->SetTrialPoints(seeds);
	//fastMarching->set
	fastMarching->SetStoppingValue(stoppingTime);
	std::cout << "Performing FastMarching Segmentation..." << std::endl;
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

	std::cout <<"Spacing		 : "<< fastMarching->GetOutput()->GetSpacing() << std::endl;
	std::cout << "SpeedConstant  : "<< fastMarching->GetSpeedConstant() << std::endl;
//	std::cout << fastMarching->GetLargeValue << std::endl;
	std::cout <<"StoppingValue	 :" << fastMarching->GetStoppingValue() << std::endl;
	typedef itk::CastImageFilter<InternalImageType, OutputImageType> castType;
	castType::Pointer castTypeFilter = castType::New();
	castTypeFilter->SetInput(fastMarching->GetOutput());
	string temp = "_FastMarching.nii.gz";
	string file = filename + temp;
	Write(castTypeFilter->GetOutput(), file);
	file.clear();  temp.clear();
	//typedef itk::MinimumMaximumImageCalculator <OutputImageType>
		//ImageCalculatorFilterType;

	//ImageCalculatorFilterType::Pointer imageCalculatorFilter
	//	= ImageCalculatorFilterType::New();
	//imageCalculatorFilter->SetImage(castTypeFilter->GetOutput());
	//imageCalculatorFilter->Compute();
	//cout << imageCalculatorFilter->GetMaximum() << endl;
	/*typedef itk::InvertIntensityImageFilter< OutputImageType, OutputImageType > InvertType;
	InvertType::Pointer inverter = InvertType::New();
	inverter->SetInput(castTypeFilter->GetOutput());
	inverter->SetMaximum(-1);
	try {
		inverter->Update();
	}
	catch (itk::ExceptionObject & excep)
	{
		std::cerr << "Exception caught in fast marching !" << std::endl;
		std::cerr << excep << std::endl;
	}

	string temp1 = "_inverstedFast.nii.gz";
	Write(inverter->GetOutput(), filename + temp1);
*/
	return fastMarching->GetOutput();
}
static OutputImageType::Pointer ThresholdFilterFast(InternalImageType::Pointer image, IndexType seed,
	const double stoppingTime)
{
	typedef itk::BinaryThresholdImageFilter<InternalImageType, OutputImageType> ThresholdType;
	ThresholdType::Pointer thresholdFilter = ThresholdType::New();
	//Watershed Segmentation Output as Threshold Output
	thresholdFilter->SetInput(image);
	thresholdFilter->SetInPlace(true);
	//Set Threshold filter limits
	float value = image->GetPixel(seed);
	//std::cout << value << std::endl;
	//thresholdFilter->SetLowerThreshold(-1000.0);
	thresholdFilter->SetUpperThreshold(stoppingTime*.99);
	//Set Output Values
	thresholdFilter->SetOutsideValue(0);
	thresholdFilter->SetInsideValue(1);
	//  std::cout<<"Extracting Seed Label..."<<std::endl;
	t.Start();
	//Apply Threshold
	thresholdFilter->Update();
	t.Stop();
	std::cout << "...threshold Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	myfile2 << "...threshold Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	// Compute some statistics
	return thresholdFilter->GetOutput();
}
static InternalImageType::Pointer ThresholdFilterFastGeo(InternalImageType::Pointer image, IndexType seed,
	const double stoppingTime)
{
	typedef itk::BinaryThresholdImageFilter<InternalImageType, InternalImageType> ThresholdType;
	ThresholdType::Pointer thresholdFilter = ThresholdType::New();
	//Watershed Segmentation Output as Threshold Output
	thresholdFilter->SetInput(image);
	thresholdFilter->SetInPlace(true);
	//Set Threshold filter limits
	float value = image->GetPixel(seed);
	std::cout << value << std::endl;
	//thresholdFilter->SetLowerThreshold(-1000.0);
	thresholdFilter->SetUpperThreshold((stoppingTime)*.90);
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