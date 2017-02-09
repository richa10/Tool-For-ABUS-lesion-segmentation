#include "itkMedianImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "utility.h"

static InternalImageType ::Pointer Smoothing(InternalImageType::Pointer image, float sigma)
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
		std::cerr << "exception in gaussian smoothing Filter " << std::endl;
		std::cerr << e << std::endl;
		return NULL;
	}
	std::cout << "...Smoothing Filter Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	return smoothingFilterZ->GetOutput();
}
static InternalImageType ::Pointer HeightFunction(InternalImageType::Pointer image)
{
	//Gradient Magnitude Filter (generates the height function)
	/*
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<InternalImageType, InternalImageType>GradientFilterType;
	GradientFilterType::Pointer gradienFilter = GradientFilterType::New();*/
	typedef itk::GradientMagnitudeImageFilter<InternalImageType, InternalImageType>GradientFilterType;
	GradientFilterType::Pointer gradienFilter = GradientFilterType::New();
	gradienFilter->SetInput(image);
	 //gradienFilter->SetSigma( 0.5 );   //sigma is  considered in millimeters
	 gradienFilter->SetUseImageSpacingOff();
	//Apply Gradient Magnitude Filter
	
	std::cout << "Applying Gradient Magnitude Filter ..." << std::endl;
	myfile2 << "Applying Gradient Magnitude Filter ..." << std::endl;
	try
	{
		t.Start();
		gradienFilter->Update();
		t.Stop();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in Gradient Magnitude Filter " << std::endl;
		std::cerr << e << std::endl;
		return NULL;
	}
	//	std::cout << "...Gradient Magnitude Filter Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	myfile2 << "...Gradient Magnitude Filter Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	return gradienFilter->GetOutput();
}
static InternalImageType::Pointer sigmoid(InternalImageType::Pointer image, float alpha, float beta, char* name)
{
	typedef   itk::SigmoidImageFilter< InternalImageType, InternalImageType >  SigmoidFilterType;
	SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
	sigmoid->SetOutputMinimum(0);
	sigmoid->SetOutputMaximum(1);
	sigmoid->SetAlpha(alpha);
	sigmoid->SetBeta(beta);
	sigmoid->SetInput(image);
	std::cout << "Applying Sigmoid Filter ..." << std::endl;
	myfile2 << "Applying Sigmoid Filter ..." << std::endl;
	try
	{
		t.Start();
		sigmoid->Update(); 
		t.Stop();
	}
	catch (itk::ExceptionObject & e) {
		std::cerr << "Sigmoid Filter Error" << std::endl;
		std::cerr << e << std::endl;
	}
	
//	std::cout << "...Sigmoid filter Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	myfile2 << "...Sigmoid filter Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	string temp1 = "_sigmoid.nii.gz";
	string file = name + temp1;
	//InternalImageType::Pointer sigmoidimage = sigmoid->GetOutput();
	Write(sigmoid->GetOutput(), file);
	return sigmoid->GetOutput();
}
static  InternalImageType::Pointer MedianSmoothing(InternalImageType::Pointer iImage, int iRadius)
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
	std::cout << "Applying Median Filter ..." << radius<< std::endl;

	try
	{
		t.Start();
		std::cout << "Applying Median Filter ..." << radius << std::endl;
		medianFilter->Update();
		std::cout << "Applying Median Filter ..." << radius << std::endl;
		t.Stop();
	}
	catch (itk::ExceptionObject & e) {
		std::cerr << "Median Filter Error" << std::endl;
		std::cerr << e << std::endl;
	}
	cout << "...median filter Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	return medianFilter->GetOutput();
}
