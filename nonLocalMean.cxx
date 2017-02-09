#include "itkPatchBasedDenoisingImageFilter.h"
#include "itkGaussianRandomSpatialNeighborSubsampler.h"
#include "../seg/utility.h"

static InternalImageType::Pointer nlmeans(InternalImageType::Pointer image
	,const int numIterations, const int numThreads, const int numToSample, const float kernelBandwidthMultiplicationFactor
	, const std::string & noiseModel,	const float noiseModelFidelityWeights )
{
	typedef   itk::PatchBasedDenoisingImageFilter
					< InternalImageType, 
						InternalImageType >  nlFilterType;
	
	typedef itk::Statistics::GaussianRandomSpatialNeighborSubsampler
					<typename nlFilterType::PatchSampleType, 
					typename InternalImageType::RegionType> SamplerType;
	
	typedef typename nlFilterType::OutputImageType OutputImageType;

	nlFilterType::Pointer nl = nlFilterType::New();
	
	nl->SetInput(image);
	// patch radius is same for all dimensions of the image
	const unsigned int patchRadius = 4;
	nl->SetPatchRadius(patchRadius);
	// instead of directly setting the weights, could also specify type
	nl->UseSmoothDiscPatchWeightsOn();
	nl->UseFastTensorComputationsOn();
	// noise model to use
	if (noiseModel == "GAUSSIAN")
	{
		nl->SetNoiseModel(nlFilterType::GAUSSIAN);
	}
	else if (noiseModel == "RICIAN")
	{
		nl->SetNoiseModel(nlFilterType::RICIAN);
	}
	else if (noiseModel == "POISSON")
	{
		nl->SetNoiseModel(nlFilterType::POISSON);
	}
	else
	{
		nl->SetNoiseModel(nlFilterType::NOMODEL);
	}
	// stepsize or weight for smoothing term
	// Large stepsizes may cause instabilities.

	nl->SetSmoothingWeight(1);

	// stepsize or weight for fidelity term
	// use a positive weight to prevent oversmoothing
	// (penalizes deviations from noisy data based on a noise model)
	nl->SetNoiseModelFidelityWeight(noiseModelFidelityWeights);

	// number of iterations over the image of denoising
	nl->SetNumberOfIterations(numIterations);

	// number of threads to use in parallel
	nl->SetNumberOfThreads(numThreads);

	// sampling the image to find similar patches
	typename SamplerType::Pointer sampler = SamplerType::New();

	// variance (in physical units) for semi-local Gaussian sampling
	//sampler->SetVariance(100);

	//// rectangular window restricting the Gaussian sampling
	//sampler->SetRadius(10); // 2.5 * standard deviation
	//						// number of random sample "patches" to use for computations
	//sampler->SetNumberOfResultsRequested(numToSample);

	//// Sampler can be complete neighborhood sampler, random neighborhood sampler,
	//// Gaussian sampler, etc.
	//nl->SetSampler(sampler);

	// automatic estimation of the kernel bandwidth
	nl->KernelBandwidthEstimationOn();

	// update bandwidth every 'n' iterations
	nl->SetKernelBandwidthUpdateFrequency(3);

	// use 33% of the pixels for the sigma update calculation
	nl->SetKernelBandwidthFractionPixelsForEstimation(0.20);

	// multiplication factor modifying the automatically-estimated kernel sigma
	nl->SetKernelBandwidthMultiplicationFactor(kernelBandwidthMultiplicationFactor);

	// manually-selected Gaussian kernel sigma
	// filter->DoKernelBandwidthEstimationOff();
	// typename FilterType::RealArrayType gaussianKernelSigma;
	// gaussianKernelSigma.SetSize(reader->GetOutput()->GetNumberOfComponentsPerPixel());
	// gaussianKernelSigma.Fill(11);
	// filter->SetGaussianKernelSigma (gaussianKernelSigma);
	try
	{
		nl->Update();
	}
	catch (itk::ExceptionObject & excp)
	{
		std::ostringstream itkmsg;                                            \
			itkmsg << "Error: In " __FILE__ ", line " << __LINE__ << "\n"
			<< "Caught exception <" << excp
			<< "> while running patch-based denoising image filter."
			<< "\n\n";
		::itk::OutputWindowDisplayWarningText(itkmsg.str().c_str());
		return ITK_NULLPTR;
	}

	
	return  nl->GetOutput();
}


//gtImageout = tempObj.ReadImageOut(argv[6]);                                                                                                                           //load the ground truth
//gtImage = tempObj.ReadImage(argv[6]);
//gtImageN = tempObj.normalizeGtImage(gtImage);

//unsigned int numIterations = 1;
//if (argc > 8)
//{
//	numIterations = atoi(argv[8]);
//}

//unsigned int numThreads = 1;
//if (argc > 9)
//{
//	numThreads = atoi(argv[9]);
//}

//unsigned int numToSample = 1000;
//if (argc > 10)
//{
//	numToSample = atoi(argv[10]);
//}

//float kernelBandwidthMultFactor = 1;
//if (argc > 11)
//{
//	kernelBandwidthMultFactor = atof(argv[11]);
//}

//std::vector< std::string > modelChoices;
//modelChoices.push_back("GAUSSIAN");
//modelChoices.push_back("RICIAN");
//modelChoices.push_back("POISSON");
//modelChoices.push_back("NOMODEL");
//std::string noiseModel;
//noiseModel = modelChoices[0];
//
//float noiseModelFidelityWeight = 0.0;
//if (argc > 12)
//{
//	noiseModel = argv[12];
//	bool validChoice = false;
//	for (unsigned int ii = 0; ii < modelChoices.size(); ++ii)
//	{
//		if (noiseModel == modelChoices[ii])
//		{
//			validChoice = true;
//		}
//	}
//	if (!validChoice)
//	{
//		std::cerr << noiseModel << " is not a valid noise model choice.  Please choose one of: ";
//		for (unsigned int ii = 0; ii < modelChoices.size(); ++ii)
//		{
//			std::cerr << modelChoices[ii] << " " << std::endl;
//		}
//		return EXIT_FAILURE;
//	}
//	if (argc > 13)
//	{
//		noiseModelFidelityWeight = atof(argv[13]);
//	}
//	else
//	{
//		std::cerr << "Must also specify a noise model fidelity weight when a noise model is specified."
//			<< std::endl;
//		return EXIT_FAILURE;
//	}
//}
//Truncating the whole image to get a mask of 70mm
/*sss << "_Mask.mha"; stn = sss.str(); newName = name + stn;
char *tempname = new char[newName.length() + 1];
std::strcpy(tempname, newName.c_str());
cout << tempname << endl;*/

/*maskImageNew = segment.CreateMask(diacomImage, diacomImage, axesImage, seedPoint);

segment.Write(maskImageNew, newName);

newName.clear(); sss.str(""); stn.clear();*/
//	string b = "rescled.mha";

//InternalImageType:: Pointer rescalImage = segment.RescaleImage(diacomImage);
//cout << "performing Non local mean" << endl;
//cout << numIterations << ", " << numThreads << ", " << numToSample << ",  " << kernelBandwidthMultFactor << ", " << noiseModelFidelityWeight << ", " << noiseModel << endl;
//nlmean.nlmeans(rescalImage, numIterations, numThreads, numToSample, kernelBandwidthMultFactor, noiseModel, noiseModelFidelityWeight);
