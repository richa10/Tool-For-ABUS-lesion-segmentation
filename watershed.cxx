#include "itkWatershedImageFilter.h"
#include "utility.h"


static  OutputImageType::Pointer watershedSegmentation(InternalImageType::Pointer image, float threshold,
	float level, IndexType seed, char *filename)
{
	//Watershed segmentation
	typedef itk::WatershedImageFilter<InternalImageType> WatershedFilterType;
	WatershedFilterType::Pointer watershedFilter = WatershedFilterType::New();
	WatershedImageType::RegionType region;
	//Gradient Magnitude Filter Output as Watershed Input
	watershedFilter->SetInput(image);
	//Threshold
	watershedFilter->SetThreshold(threshold);
	//Level
	watershedFilter->SetLevel(level);
	std::cout << "Performing Watershed Segmentation" << std::endl;
	myfile2 << "Performing Watershed Segmentation" << std::endl;
	myfile2 << "seedPoint in watershed  " << seed[0] << "," << seed[1] << "," << seed[2] << "\n";

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
	//std::cout << "the tresholded intensity is " << value << endl;
	myfile2 << "the tresholded intensity is " << value << endl;
	float newlevel = watershedFilter->GetLevel();;
	if (value == 1)
	{
		while (value <= 1 && newlevel > 0)
		{
			newlevel = newlevel - 0.1;
			myfile2 << "new level value" << newlevel << endl;
			myfile2 << "Performing Watershed Segmentation." << std::endl;
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
			myfile2 << "the tresholded intensity is " << value << endl;
			newlevel = watershedFilter->GetLevel();
		}
		newlevel = newlevel + 0.1;
		value = 1;
		while (value <= 1 && newlevel > 0)
		{
			newlevel = newlevel - 0.05;
			myfile2 << "new level value" << newlevel << endl;
			myfile2 << "Performing Watershed Segmentation.." << std::endl;
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
			myfile2 << "the tresholded intensity is " << value << endl;
			newlevel = watershedFilter->GetLevel();
		}
		
		newlevel = newlevel + 0.05;
		value = 1;
		while (value <= 1 && newlevel > 0)
		{
			newlevel = newlevel - 0.01;
			myfile2 << "new level value" << newlevel << endl;
			myfile2 << "Performing Watershed Segmentation..." << std::endl;
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
			myfile2 << "the tresholded intensity is " << value << endl;
			newlevel = watershedFilter->GetLevel();
		}
		myfile2 << "the tresholded intensity is " << value <<"and level "<< newlevel << endl;
		//cout << "the tresholded intensity is " << value << endl;
	}
	//std::cout << "...Watershed segmentation Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	myfile2 << "...Watershed segmentation Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;

	std::stringstream ss;
	ss << "_" << newlevel << "_" << value << "_" <<seed[0]+1<<"_" << seed[1]+1 << "_"<< seed[2]+1 << "_";
	std::string temp = ss.str();
	string waterShedfile = filename + temp + "watershed.nii.gz";
	Write(watershedFilter->GetOutput(),waterShedfile);
	//region = watershedFilter->GetOutput()->get
	//typedef itk::CastImageFilter< InternalImageType, OutputImageType > CastingFilterType2;
	//CastingFilterType2::Pointer Castingfilter2 = CastingFilterType2::New();
	//Castingfilter2->SetInput(watershedFilter->GetOutput());
	//Castingfilter2->Update();
	////Seed Label
	//// WriteDiacomImage(watershedFilter->GetOutput(),waterShedfile);
	//try
	//{
	//	toolsITK::writeDCMITKImage(Castingfilter2->GetOutput(), waterShedfile);
	//}
	//catch (itk::ExceptionObject & e)
	//{
	//	std::cerr << "exception in file writer " << std::endl;
	//	std::cerr << e << std::endl;
	//	return 0;
	//}
	OutputPixelType valueInten = watershedFilter->GetOutput()->GetPixel(seed);
	//cout << "the label value is   " << valueInten << endl;
	myfile2 << "the label value is   " << valueInten << endl;
	//Label Threshold
	typedef itk::BinaryThresholdImageFilter<WatershedImageType, OutputImageType> ThresholdType;
	ThresholdType::Pointer thresholdFilter = ThresholdType::New();
	//Watershed Segmentation Output as Threshold Output
	thresholdFilter->SetInput(watershedFilter->GetOutput());
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
	//std::cout << "...threshold Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;
	myfile2 << "...threshold Success(" << t.GetMean() << " s./" << t.GetMean() / 60.0 << " min.)" << std::endl;

	//// Compute some statistics
	//double sourceImageSpacing[] = { image->GetSpacing()[0],image->GetSpacing()[1],image->GetSpacing()[2] };  // (pixels)
	//																										 //std::cout<< "Image spacing: " << sourceImageSpacing[0] << " x " << sourceImageSpacing[1] << " x " << sourceImageSpacing[2] << " mm^3" << std::endl;
	//getVolume(thresholdFilter->GetOutput(), sourceImageSpacing);
	//ss.str(""); temp.clear(); waterShedfile.clear();
	//ss << "_" << iterations; temp = ss.str();
	//waterShedfile = filename + temp + "_segmented.dcm";
	//try
	//{
	//	toolsITK::writeDCMITKImage(thresholdFilter->GetOutput(), waterShedfile);
	//}
	//catch (itk::ExceptionObject & e)
	//{
	//	std::cerr << "exception in file writer " << std::endl;
	//	std::cerr << e << std::endl;
	//	return NULL;
	//}
	return thresholdFilter->GetOutput();
}