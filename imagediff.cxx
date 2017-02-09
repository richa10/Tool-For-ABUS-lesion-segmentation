#include "itkImage.h"
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


typedef itk::Image< unsigned char, 3>  OutImageType;
typedef itk::Image<float, 3>  FloatImageType;

int main(int argc, char *argv[])
{
	typedef itk::ImageFileReader<FloatImageType> ReaderType;  //the diacom image
	ReaderType::Pointer readImage1 = ReaderType :: New();
	ReaderType::Pointer readImage2 = ReaderType::New();
	readImage1->SetFileName(argv[1]);
	readImage2->SetFileName(argv[2]);
	
	try {
		readImage1->Update();
		readImage2->Update();
	}
	catch (itk::ExceptionObject & e) {
		std::cerr << "reader Filter Error" << std::endl;
		std::cerr << e << std::endl;
	}
	typedef itk::AbsoluteValueDifferenceImageFilter <FloatImageType, FloatImageType,
		OutImageType>
          AbsoluteValueDifferenceImageFilterType;
  AbsoluteValueDifferenceImageFilterType::Pointer absoluteValueDifferenceFilter
          = AbsoluteValueDifferenceImageFilterType::New ();
  absoluteValueDifferenceFilter->SetInput1(readImage1->GetOutput());
  absoluteValueDifferenceFilter->SetInput2(readImage2->GetOutput());
  
  try {
	
	  absoluteValueDifferenceFilter->Update();
  }
  catch (itk::ExceptionObject & e) {
	  std::cerr << "difference Filter Error" << std::endl;
	  std::cerr << e << std::endl;
  }

  typedef itk::ImageFileWriter<OutImageType> WriterType;
  WriterType::Pointer writeImage1 = WriterType::New();
  writeImage1->SetFileName("G:/codeAnisotropic/result_.06_nl/FileDifference1.mha");
  writeImage1->SetInput(absoluteValueDifferenceFilter->GetOutput());
  
  try {
	    writeImage1->Update();
  }
  catch (itk::ExceptionObject & e) {
	  std::cerr << "writer Filter Error" << std::endl;
	  std::cerr << e << std::endl;
  }
 
  return EXIT_SUCCESS;
}


