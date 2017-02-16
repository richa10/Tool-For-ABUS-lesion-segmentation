
//#include "../seg/nonLocalMean.cxx"
#include "AnisotropicLBR.cxx"
#include "geodesic.cxx"
#include "watershed.cxx"
#include "createMask.cxx"
#include "fastmarching.cxx"
#include "preprocessing.cxx"
#include "math.h"
#include "utility.h"
#include "itkImageMomentsCalculator.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
//#include "itkInvertIntensityImageFilter.h.";
int main(int argc, char *argv[])
{
	//cout << argc << endl;
	if (argc != 16) {
		cout << "Input_image Output_image SeedX SeedY SeedZ" << std::endl;
		return 0;
	}
	itk::TimeProbe clock;
	clock.Start();
	char * name = argv[5];
	char * file = argv[7];
	std::string newName;
	float *seedPointPhy = 0;
	seedPointPhy = new float[3];
	seedPointPhy[0] = atof(argv[2])+ atof(argv[15]); // (pixels)
	seedPointPhy[1] = atof(argv[3])- atof(argv[15]); // (pixels)
	seedPointPhy[2] = atof(argv[4])- atof(argv[15]); // (pixels)

	IndexType axesImage;
	axesImage[0] = 75; // atoi(argv[16])*2;      // diameter (mm)
	axesImage[1] = 75;// atoi(argv[17])*2;      // (mm)
	axesImage[2] = 75; //atoi(argv[18])*2;      // (mm)

	IndexType axes;
	axes[0] = axesImage[0] / 2;      // (mm)
	axes[1] = axesImage[0] / 2;      // (mm)
	axes[2] = axesImage[0] / 2;      // (mm)

	IndexType axesSeed;
	axesSeed[0] = 5;      // diameter (mm)
	axesSeed[1] = 5;      // (mm)
	axesSeed[2] = 5;      // (mm)

	float* sigMuStd = 0;
	sigMuStd = new float[2];
	float* sigMuStdseed = 0;
	sigMuStdseed = new float[2];

	//to store dice and other evaluation measure
	float *labelMeasure = 0;
	labelMeasure = new float[3];

	//Perform algorithm
	InternalImageType::Pointer  diacomImage, gtImage, smoothimage;
	OutputImageType::Pointer segmentationImageGeo, segmentationImageFastMarch, gtImageout, segmentationImageWater;
	std::stringstream sss;
	std::string stn;



	myfile2.open(file, ios::out | ios::app);
	diacomImage = ReadImage(argv[1]);
	double sourceImageSpacing[] = { diacomImage->GetSpacing()[0],diacomImage->GetSpacing()[1],diacomImage->GetSpacing()[2] }; //get the spacing of the image
	double origin[] = { diacomImage->GetOrigin()[0],diacomImage->GetOrigin()[1], diacomImage->GetOrigin()[2] };
	IndexType seedPoint;
	seedPoint[0] = floor(((seedPointPhy[0] - origin[0]) / sourceImageSpacing[0]) + 0.5);
	seedPoint[1] = floor(((seedPointPhy[1] - origin[1]) / sourceImageSpacing[1]) + 0.5);
	seedPoint[2] = floor(((seedPointPhy[2] - origin[2]) / sourceImageSpacing[2]) + 0.5);
	


	/*myfile2 << "seedPoint  " << seedPoint[0] << "," << seedPoint[1] << "," << seedPoint[2] << "\n";
	cout << "seedPoint  " << seedPoint[0] << "," << seedPoint[1] << "," << seedPoint[2] << "\n";*/

	gtImage = ReadImage(argv[6]);

	string sb = "";
	gtImageout = ThresholdFilterGT(gtImage);
	
	InternalImageType::Pointer MaskImage = CreateMask(diacomImage, diacomImage, axesImage, seedPoint);
	sigMuStdseed = getMuSigma(MaskImage);
	
	int	diffusionTime = atoi(argv[8]);
	double  lambda = atof(argv[9]);
	char* enhancement = argv[10];
	double	featureScale = atoi(argv[11]);
	double exponent = atoi(argv[12]);
	double noiseScale = 0.6;// sigMuStdseed[1] * 0.025;
	double diffparam = (sigMuStdseed[0] - sigMuStdseed[1])*0.01;
	double snr = sigMuStdseed[0] / sigMuStdseed[1];
	
	diffusionTime = 1.5; //0.2 mm images
	if (snr < 1.5)
	{
		noiseScale = 0.3;
		lambda = 0.03;
		diffusionTime = 1.0;
		enhancement = "cCED";//cCED
		//featureScale = 4; //4
	}
	else if (1.5 <= snr && snr < 1.8)
	{
		noiseScale = 0.8;//1.0
		lambda = 0.003;  //0.001
		diffusionTime = 2; //1.5
		enhancement = "cEED"; //CED
		//featureScale = 3; //4
	}
	else if (1.8 <= snr && snr < 2)
	{
		noiseScale = 0.7;//1.0
		lambda = 0.001;  //0.001
		diffusionTime = 1.5; //1.5
		enhancement = "cEED"; //CED
							  //featureScale = 3; //4
	}
	else if (2 <= snr && snr < 2.4)
	{

		noiseScale = 0.8; //0.5
		lambda = 0.005;  //0.005
		diffusionTime =1.5; //1.5
		enhancement = "cEED"; //cEED
		//featureScale = 3; //4
		cout << enhancement;
	}
	else if (2.4 <= snr && snr < 2.6)
	{
		noiseScale = 0.8; //0.8
		lambda = 0.005;  //0.005
		diffusionTime = 1.5; //1.5
		enhancement = "CED"; //CED
							 //featureScale = 4; //4
		cout << enhancement;
	}
	else if (2.6 <= snr && snr <3.1)
	{
		noiseScale = 0.6; //0.8
		lambda = 0.005;  //0.005
		diffusionTime = 1.5; //1.5
		enhancement = "cEED"; //CED
		featureScale = 3; //4
		cout << enhancement;
	}
	
	else if (3.1 <= snr && snr < 3.4)
	{
		noiseScale = 0.6; //0.6
		lambda = 0.005;  //0.005
		diffusionTime = 2.0; //1.5
		enhancement = "cEED"; //cEED
		//featureScale = 4; //4
		cout << enhancement;
	}
	else if (3.4 <= snr && snr < 4.0)
	{
		noiseScale = 0.6;//1//0.6 //0.2
		lambda = 0.5;//0.5//0.05 //0.005
		diffusionTime = 1.0;//2 //1
		enhancement = "cEED"; //cCED
		//featureScale = 3;
		//exponent = 3;
	}
	else
	{
		noiseScale = 0.5;//1//0.6
		lambda = 0.0005;//0.5//0.05//0.00005
		diffusionTime = 1.5;//2
		enhancement = "CED";//cEED
		//featureScale = 3;
		//exponent = 3;

	}
	
	smoothimage = coherenceAnisotropic(MaskImage, diffusionTime, lambda, enhancement, noiseScale, featureScale, exponent);
	myfile2 << "stats of non local image output" << std::endl;
	string s = enhancement;
	stn = "_.nii.gz";
	Write(MaskImage, name + s + stn);
	stn.clear();
	MaskImage = NULL;
	sigMuStd[0] = 0; sigMuStd[1] = 0;
	MaskImage = CreateMask(smoothimage, smoothimage, axesImage, seedPoint);
	sigMuStd = getMuSigma(MaskImage);
	myfile2 << "after smoothing ###############################------------" << ", " << ", " << sigMuStd[0] << " , " << sigMuStd[1] << std::endl;
	//cout << "after smoothing ###############################------------" << ", " << ", " << sigMuStd[0] << " , " << sigMuStd[1] << std::endl;


	//calculating the image and gt gradient image or height function
	InternalImageType::Pointer  heightImage;
	//stn = "_GM.nii.gz";
	heightImage = HeightFunction(smoothimage);
	//Write(heightImage, name + stn);
	stn.clear();

	// creating the mask again for the image to avoide the edges artifact due to masking
	InternalImageType::Pointer maskSmoothImage;
	//stn = "_GM_Mask.nii.gz";       newName = name + stn;
	maskSmoothImage = CreateMask(diacomImage, heightImage, axes, seedPoint);
	//Write(maskSmoothImage, newName);
	newName.clear(); stn.clear();

	//Water shed segmentation
	segmentationImageWater = watershedSegmentation(maskSmoothImage, 0.005, 0.8, seedPoint, name);
	string wa = "_watershedSgmentation.nii.gz";
	Write(segmentationImageWater, name + wa);
	float volumeGt = getVolume(gtImageout, sourceImageSpacing);

	float volumeSegWater = getVolume(segmentationImageWater, sourceImageSpacing);

	float* vecWater = 0;
	vecWater = new float[3];
	vecWater = evaluate(gtImageout, segmentationImageWater);
	myfile2 << "Union (jaccard)" << "," << vecWater[0] << "\n";
	myfile2 << "Mean (dice)" << "," << vecWater[1] << "\n";
	myfile2 << "False negative" << "," << vecWater[2] << "\n";
	myfile2 << "False positive" << "," << vecWater[3] << "\n";

	clock.Stop();

	myfile2.close();

	myfile2.open(argv[14], ios::out | ios::app);
	myfile2 << argv[13] << "," << volumeGt << "," << snr << "," << ",";
	myfile2 << vecWater[0] << "," << vecWater[1] << "," << vecWater[2] << "," << vecWater[3] << "," << volumeSegWater << "," << endl;
	
	myfile2.close();


	// print computation time at screen
	myfile2.open(file, ios::out | ios::app);
	myfile2 << "###############################" << std::endl;
	myfile2 << "SegmentationABUS" << std::endl;
	myfile2 << "Total time: " << clock.GetTotal() / 60.0 << " min." << std::endl;
	myfile2 << "Total time: " << clock.GetTotal() << " s." << std::endl;
	myfile2 << "###############################" << std::endl;
	myfile2 << "***************************************************************" << std::endl << endl;
	myfile2.close();

	cout << "###############################" << std::endl;
	cout << "SegmentationABUS" << std::endl;
	cout << "Total time: " << clock.GetTotal() / 60.0 << " min." << std::endl;
	cout << "Total time: " << clock.GetTotal() << " s." << std::endl;
	cout << "###############################" << std::endl;
	cout << "***************************************************************" << std::endl << endl;
}
