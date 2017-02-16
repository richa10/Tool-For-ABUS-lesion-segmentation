
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkAnisotropicDiffusionImageFilter.h"
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include<itkCoherenceEnhancingDiffusionFilter.h>
#include "utility.h"

//coherenceAnisotropic(diacomImage, numIterations, numThreads, numToSample, kernelBandwidthMultFactor, noiseModel, noiseModelFidelityWeight);
static  InternalImageType::Pointer coherenceAnisotropic(InternalImageType::Pointer image, const float diffusionTime, const double lambda,
	const char * enhancement, const double noiseScale, const double featureScale, const double exponent)
{

	//typedef itk::AnisotropicDiffusionLBRImageFilter<InternalImageType, float>
	typedef itk::CoherenceEnhancingDiffusionFilter<InternalImageType, float> DiffusionFilterType;
	DiffusionFilterType::Pointer diffusionFilter = DiffusionFilterType::New();
	diffusionFilter->SetInput(image);

	diffusionFilter->SetDiffusionTime(diffusionTime);
	diffusionFilter->SetLambda(lambda);
	diffusionFilter->SetAlpha(0.01);
	
	if (!strcmp(enhancement, "EED"))
		diffusionFilter->SetEnhancement(DiffusionFilterType::EED); // Weickert's exponent : 4.
	else if (!strcmp(enhancement, "cEED"))
		diffusionFilter->SetEnhancement(DiffusionFilterType::cEED); // Weickert's exponent : 4.
	else if (!strcmp(enhancement, "CED"))
		diffusionFilter->SetEnhancement(DiffusionFilterType::CED); // Weickert's exponent : 2.
	else if (!strcmp(enhancement, "cCED"))
		diffusionFilter->SetEnhancement(DiffusionFilterType::cCED); // Weickert's exponent : 2.
	else if (!strcmp(enhancement, "Isotropic"))
		diffusionFilter->SetEnhancement(DiffusionFilterType::Isotropic); //Perona-Mali's exponent: 2.


	

	diffusionFilter->SetNoiseScale(noiseScale);
	diffusionFilter->SetFeatureScale(featureScale);
	
	
	diffusionFilter->SetExponent(exponent);
	diffusionFilter->SetAdimensionize(true);
	//diffusionFilter->
	
	/*cout <<
	"T: " << diffusionFilter->GetDiffusionTime() << "\n" <<
	"Lambda: " << diffusionFilter->GetLambda() << "\n"  <<
	"enhancement: " << diffusionFilter->GetEnhancement() << "  "<< enhancement <<"\n";
*/
	
	try
	{
		t.Start();
		diffusionFilter->Update();
		t.Stop();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception caught here in lbr!! !" << std::endl;
		std::cerr << e << std::endl;
		return NULL;
	}
	//std::cout << "...AnisotropicLBR diffusion Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	myfile2 << "...AnisotropicLBR diffusion Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;

	auto imageout = diffusionFilter->GetOutput();

	
	myfile2<< "GetLinearFilterEffectiveTimesAndIterations  "<< diffusionFilter->GetLinearFilterEffectiveTimesAndIterations()[0].first <<"  " << diffusionFilter->GetLinearFilterEffectiveTimesAndIterations()[1].second << std::endl;
	return diffusionFilter->GetOutput();

	}

static InternalImageType::Pointer anisotropicDiffusion(InternalImageType::Pointer imagei, float conductance, int iterations)
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
	//must be smaller than  = imagespacing / 2^(dimention +1) used is step  = imagespacing / 2^(dimention +1)
	float step = imagei->GetSpacing()[0] / pow(2, Dimension + 1);
	step = step*.99;
	//std::cout << "the step is :: " << step <<"\n";
	anisotropic->SetTimeStep(step);
	//Conductance
	anisotropic->SetConductanceParameter(conductance);
	//std::cout << "Applying Anistropic Diffusion ..." << std::endl;
	myfile2 << "Applying Anistropic Diffusion ..." << std::endl;
	//Apply Anisotropic Difussion
	t.Start();
	anisotropic->Update();
	t.Stop();
	//std::cout << "...Anisotropic diffusion Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	myfile2 << "...Anisotropic diffusion Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	return anisotropic->GetOutput();
}
