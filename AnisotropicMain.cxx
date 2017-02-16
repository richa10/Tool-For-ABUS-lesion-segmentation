//#include "../seg/nonLocalMean.cxx"
#include "AnisotropicLBR.cxx"
#include "geodesic.cxx"
#include "watershed.cxx"
#include "createMask.cxx"
#include "fastmarching.cxx"
#include "preprocessing.cxx"
#include "math.h"
#include "utility.h"

#include "itkSignedDanielssonDistanceMapImageFilter.h"
int main(int argc, char *argv[])
{
	//cout << argc << endl;
	if (argc != 10) {
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
	seedPointPhy[0] = atof(argv[2]); // (pixels)
	seedPointPhy[1] = atof(argv[3]); // (pixels)
	seedPointPhy[2] = atof(argv[4]); // (pixels)

	IndexType axesImage;
	axesImage[0] = 75; // atoi(argv[16])*2;      // diameter (mm)
	axesImage[1] = 75;// atoi(argv[17])*2;      // (mm)
	axesImage[2] = 75; //atoi(argv[18])*2;      // (mm)

	IndexType axes;
	axes[0] = axesImage[0] / 2;      // (mm)
	axes[1] = axesImage[0] / 2;      // (mm)
	axes[2] = axesImage[0] / 2;      // (mm)

	IndexType axesSeed;
	axesSeed[0] = 7;      // diameter (mm)
	axesSeed[1] = 7;      // (mm)
	axesSeed[2] = 7;      // (mm)

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
	myfile2 << "\n" << argv[5] << "\n" << endl;
	myfile2 << "physical co-ordinates " << atof(argv[2]) << "," << atof(argv[3]) << "," << atof(argv[4]) << "\n";
	myfile2 << "spacing " << sourceImageSpacing[0]<<", " << sourceImageSpacing[1] << ", " << sourceImageSpacing[2] << "\n";
	myfile2 << "Origin" << origin[0] << "," << origin[1] << "," << origin[2]  << "\n";
	cout << "\n" << argv[5] << "\n" << endl;
	cout << "physical co-ordinates " << atof(argv[2]) << "," << atof(argv[3]) << "," << atof(argv[4]) << "\n";
	cout << "spacing " << sourceImageSpacing[0] << ", " << sourceImageSpacing[1] << ", " << sourceImageSpacing[2] << "\n";
	cout << "Origin" << origin[0] << "," << origin[1] << "," << origin[2] << "\n";

	IndexType seedPoint;
	seedPoint[0] = floor(((seedPointPhy[0] - origin[0]) / sourceImageSpacing[0]) + 0.5);
	seedPoint[1] = floor(((seedPointPhy[1] - origin[1]) / sourceImageSpacing[1]) + 0.5);
	seedPoint[2] = floor(((seedPointPhy[2] - origin[2]) / sourceImageSpacing[2]) + 0.5);
	
	
	myfile2 << "seedPoint  " << seedPoint[0] << "," << seedPoint[1] << "," << seedPoint[2] << "\n";
	cout << "seedPoint  " << seedPoint[0] << "," << seedPoint[1] << "," << seedPoint[2] << "\n";
	
	gtImage = ReadImage(argv[6]);
	myfile2 << "Intensity at seed in GT before normalization  " << gtImage->GetPixel(seedPoint) << endl;
	cout << "Intensity at seed in GT before normalization  " << gtImage->GetPixel(seedPoint) << endl;
	string sb =""; 
	gtImageout = ThresholdFilterGT(gtImage);
	sb = "_GT_threshold.nii.gz";
	Write(gtImageout, name + sb);
	cout << "Intensity at seed in GT " << gtImageout->GetPixel(seedPoint) << endl;
	myfile2 << "Intensity at seed in GT " << gtImageout->GetPixel(seedPoint) << endl;
	InternalImageType::Pointer MaskImage = CreateMask(diacomImage, diacomImage, axesImage, seedPoint);
	stn = "mask.nii.gz";
	Write(MaskImage, name + stn);
	sigMuStdseed = getMuSigma(MaskImage);
	myfile2 << "75 mm mask's mean and std ###############################------------"<< "," <<"," << sigMuStdseed[0] << " , " << sigMuStdseed[1] << std::endl;
	cout << "75 mm  mask's mean and std ###############################------------" << "," << "," << sigMuStdseed[0] << " , " << sigMuStdseed[1] << std::endl;
	
	float conductance =3;
	int iterations =20;

	cout << "conductance::" << conductance << endl;
	cout << "iterations :: " << iterations << endl;
	myfile2 << "conductance::" << conductance << endl;
	myfile2 << "iterations :: " << iterations << endl;
	smoothimage = anisotropicDiffusion(MaskImage, conductance, iterations);
	
	stn = "_AD.nii.gz";       newName = name + stn;
	Write(smoothimage, newName);
	newName.clear(); stn.clear();

	MaskImage = NULL;
	sigMuStd[0] = 0; sigMuStd[1] = 0;
	MaskImage = CreateMask(smoothimage, smoothimage, axesImage, seedPoint);
		
	//calculating the image and gt gradient image or height function
	InternalImageType::Pointer  heightImage;
	stn = "_GM.nii.gz";
	heightImage = HeightFunction(smoothimage);
	Write(heightImage, name + stn);
	stn.clear();

	// creating the mask again for the image to avoide the edges artifact due to masking
	InternalImageType::Pointer maskSmoothImage;
	stn = "_GM_Mask.nii.gz";       newName = name + stn;
	maskSmoothImage = CreateMask(diacomImage, heightImage, axes, seedPoint);
	Write(maskSmoothImage, newName);
	newName.clear(); stn.clear();

	//Water shed segmentation
	segmentationImageWater = watershedSegmentation(maskSmoothImage, 0.005, 0.8, seedPoint, name);
	string wa = "_watershedSgmentation.nii.gz";
	Write(segmentationImageWater, name + wa);
	cout << "Volume of the Gound truth... ";
	myfile2 << "Volume of the Gound truth...";
	float volumeGt = getVolume(gtImageout, sourceImageSpacing);

	cout << "Volume of the watershed segmentation... ";
	myfile2 << "Volume of the watershed segmentation... ";
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

	myfile2.open(argv[9], ios::out | ios::app);
	myfile2 << argv[8] << "," << volumeGt << "," << sigMuStdseed[0] << "," << sigMuStdseed[1] << "," << conductance <<", "<<iterations<< "," << ",";
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
