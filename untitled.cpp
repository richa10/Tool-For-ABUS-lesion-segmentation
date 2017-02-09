    #if defined(_MSC_VER)
    #pragma warning ( disable : 4786 )
    #endif

    #ifdef __BORLANDC__
    #define ITK_LEAN_AND_MEAN
    #endif

    #include "itkImage.h"

    #include <itkIndex.h>
    #include "itkImageFileReader.h"
    #include "itkImageFileWriter.h"
    #include "itkCastImageFilter.h"
    #include <itkNiftiImageIO.h>
    #include "itkRescaleIntensityImageFilter.h"
    #include "itkBinaryThresholdImageFilter.h"
    #include "itkStatisticsImageFilter.h"


    #include "itkGradientAnisotropicDiffusionImageFilter.h"
    #include "itkGradientMagnitudeImageFilter.h"
    #include "itkRecursiveGaussianImageFilter.h"
    #include "itkWatershedImageFilter.h"
    #include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"

    #include "itkEllipsoidInteriorExteriorSpatialFunction.h"
    #include "itkFloodFilledSpatialFunctionConditionalIterator.h"
    #include "itkMaskImageFilter.h"


    #include <stdio.h>
    #include <string>
    #include <sstream>
    #include <iostream>
    #include <math.h>
    #include <fstream>

    #include "../albertLib/toolsITK.h"
    #include "../albertLib/toolsITK.cpp"



    #include "itkTimeProbe.h"
    //Template for the Dimension
    //Image Properties
    const unsigned int  Dimension = 3;
    typedef float InternalPixelType;
    typedef itk::Image<InternalPixelType, Dimension> InternalImageType;
    typedef unsigned short OutputPixelType;
    //typedef unsigned char OutputPixelType;
    typedef itk::Image<OutputPixelType, Dimension> OutputImageType;
    typedef  typename InternalImageType::IndexType  IndexType;
    typedef  typename InternalImageType::IndexType  Indexarray[];
    typedef long unsigned int WatershedPixelType;
    typedef itk::Image<WatershedPixelType, Dimension> WatershedImageType;

    class Segmentation
    {
        itk::TimeProbe t;
        typedef itk::ImageFileReader<InternalImageType> ReaderType;  //the diacom image
        typename ReaderType::Pointer readImage;

           
        typedef itk::ImageFileWriter<WatershedImageType> WriterType;
        typename WriterType::Pointer writer;
         public: typename InternalImageType::Pointer image;
         typename InternalImageType::RegionType region;
         typename InternalImageType::SizeType size;


        //constructor takes input and output path of image

    public: Segmentation()
        {

            //instantiate reader
            readImage = ReaderType::New();

            //instantiate Writer
            writer = WriterType::New();

            //interWriter = IntermediateWriterType :: New();

        }
        //constructor takes input path of image to write

        //Read function should return InputImageType
    public: InternalImageType::Pointer ReadImage(char *filename)
        {

            readImage->SetFileName( filename);
            try
            {
                readImage->Update();
            }
            catch (itk::ExceptionObject & e)
            {
                std::cerr << "exception in file reader " << std::endl;
                std::cerr << e << std::endl;
                //return NULL;
            }
            // display input image size
            image  = readImage->GetOutput();
            region = image->GetLargestPossibleRegion();
            size = region.GetSize();
            std::cout << "Input image size: " << size << std::endl;
            return readImage->GetOutput();
        }

    public: float** getStatatistics(InternalImageType::Pointer imageStats)
        {
            // Compute some statistics
            typedef itk::StatisticsImageFilter<InternalImageType> StatisticsImageFilterType;
            typename StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
            statisticsImageFilter->SetInput(imageStats);
            statisticsImageFilter->Update();
                   double mean  = statisticsImageFilter->GetMean();
                   double sigma = statisticsImageFilter->GetSigma();
            //        double snr = (mean / sigma)*100;
            // Print statistics
            std::cout << std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
            std::cout << "Total number of voxels within lesion: " << statisticsImageFilter->GetSum() << std::endl;
            // std::cout << "Total volume of lesion: " << (statisticsImageFilter->GetSum()) <<  std::endl;
            //  std::cout << "Total Mean of lesion: " << (statisticsImageFilter->GetMean()) <<   std::endl;
            // std::cout << "Total Sigma of lesion: " << (statisticsImageFilter->GetSigma()) <<  std::endl;
            //std::cout << "SNR of the image is " << snr <<  std::endl;
            // std::cout << "Maximum of lesion: " << (statisticsImageFilter->GetMaximum()) <<  std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
            float** vec = 0;
            vec = new float*[1];
           //std::cout<<sizeof(cubeCordinates)<<endl;
           for (int h = 0; h < 1; h++)
           {
               vec[h] = new float[2];
           }
           vec[0][0] = mean;
           vec[0][1] = sigma;
           return vec;
            
        }
    public: int getStatatistics(OutputImageType::Pointer imageStats, double sourceImageSpacing[])
        {
            // Compute some statistics
            typedef itk::StatisticsImageFilter<OutputImageType> StatisticsImageFilterType;
            typename StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
            statisticsImageFilter->SetInput(imageStats);
            statisticsImageFilter->Update();

            // Print statistics
            std::cout << std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
            //std::cout << "Total number of voxels within lesion: " << statisticsImageFilter->GetSum() << std::endl;
            std::cout << "Total volume of lesion: " << (statisticsImageFilter->GetSum())*sourceImageSpacing[0]*sourceImageSpacing[1]*sourceImageSpacing[2] <<  " mm³" << std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
            return 1;
        }
    public: InternalImageType:: Pointer anisotropicDiffusion( InternalImageType :: Pointer imagei , int iterations, float conduct)
        {
            //Anisotropic Difussion (smooth the image)
            typedef itk::GradientAnisotropicDiffusionImageFilter<InternalImageType, InternalImageType> FilterType;
            typename FilterType::Pointer anisotropic = FilterType::New();

            //Input Image as Anisotropic Diffusion Input
            anisotropic->SetInput(imagei);
            anisotropic->SetInPlace(true);

            //Number of Iterations
            anisotropic->SetNumberOfIterations(iterations);
            //Step Time
            //must be smaller than  = imagespacing / 2^(dimention +1) used is step  = imagespacing / 2^(dimention +2)
            float step = image->GetSpacing()[0]/ pow(2,Dimension+2);
            //std::cout<< "the step is :: " << step ;
            anisotropic->SetTimeStep(step);
            //Conductance
            anisotropic->SetConductanceParameter(conduct);
            //	std::cout<<"Applying Anistropic Diffusion ..."<<std::endl;
            //Apply Anisotropic Difussion
            t.Start();
            anisotropic->Update();
            t.Stop();
            std::cout<<"...Anisotropic diffusion Success(" << t.GetMean() << " s./" << t.GetMean()/60.0  << " min.)" <<std::endl;

            return Rescaler(anisotropic->GetOutput());
        }
    public: InternalImageType:: Pointer Rescaler(InternalImageType :: Pointer image)
        {
            typedef itk::RescaleIntensityImageFilter<InternalImageType, InternalImageType> RescaleFilterType;
            typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
            //Image GrayValue Limits


            rescaler->SetOutputMinimum(0);
            rescaler->SetOutputMaximum(255);
            //Anisotropic Diffusion Output as Rescaler Input
            rescaler->SetInput(image);
            rescaler->SetInPlace(true);
            //	std::cout<<"Rescaling Image..."<<std::endl;
            //Rescale the Image
            t.Start();
            rescaler->Update();
            t.Stop();
            std::cout<<"...Rescale Success(" << t.GetMean() << " s./" << t.GetMean()/60.0  << " min.)" <<std::endl;

            return rescaler->GetOutput();
        }
    public: InternalImageType:: Pointer Smoothing(InternalImageType :: Pointer image, float sigma)
        {

            //Gradient Magnitude Filter (generates the height function)
            typedef itk::RecursiveGaussianImageFilter<InternalImageType,InternalImageType>smoothingFilterType;
            typename smoothingFilterType::Pointer smoothingFilterX=smoothingFilterType::New();
            typename smoothingFilterType::Pointer smoothingFilterY=smoothingFilterType::New();
            typename smoothingFilterType::Pointer smoothingFilterZ=smoothingFilterType::New();

            smoothingFilterX->SetDirection( 0 );    //0 --> X direction
            smoothingFilterY->SetDirection( 1 );    //1 --> Y direction
            smoothingFilterZ->SetDirection( 2 );    //2 --> Z direction
            smoothingFilterX->SetOrder( smoothingFilterType::ZeroOrder );
            smoothingFilterY->SetOrder( smoothingFilterType::ZeroOrder );
             smoothingFilterZ->SetOrder( smoothingFilterType::ZeroOrder );

              //Rescaler Output as  Filter Input
             smoothingFilterX->SetNormalizeAcrossScale( true );
             smoothingFilterY->SetNormalizeAcrossScale( true );
             smoothingFilterZ->SetNormalizeAcrossScale( true );

             smoothingFilterX->SetInput( image);
             smoothingFilterY->SetInput( smoothingFilterX->GetOutput() );
             smoothingFilterZ->SetInput( smoothingFilterY->GetOutput() );


            //Sigma (standard deviation of Gaussian)
            smoothingFilterX->SetSigma(sigma);  ///parameter to optimize
            smoothingFilterY->SetSigma(sigma);  ///parameter to optimize
            smoothingFilterZ->SetSigma(sigma);  ///parameter to optimize

            std::cout<<"Applying Gradient Magnitude Filter..."<< sigma<<std::endl;
            //Apply Gradient Magnitude Filter
            t.Start();
            smoothingFilterZ->Update();
            t.Stop();
            std::cout<<"...Smoothing Filter Success(" << t.GetMean() << " s./" << t.GetMean()/60.0  << " min.)" <<std::endl;
            return smoothingFilterZ->GetOutput();
        }
    public: InternalImageType:: Pointer HeightFunction(InternalImageType :: Pointer image, float sigma)
        {

            //Gradient Magnitude Filter (generates the height function)
            typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<InternalImageType,InternalImageType>GradientFilterType;
            typename GradientFilterType::Pointer gradienFilter=GradientFilterType::New();

            //Rescaler Output as Magnitude Gradient Filter Input
            gradienFilter->SetInput(image);
            gradienFilter->SetSigma( sigma ); 
            //Apply Gradient Magnitude Filter
            t.Start();
            gradienFilter->Update();
            t.Stop();
            std::cout<<"...Gradient Magnitude Filter Success(" << t.GetMean() << " s./" << t.GetMean()/60.0  << " min.)" <<std::endl;
            return gradienFilter->GetOutput();
        }

    public: InternalImageType:: Pointer CreateMask(InternalImageType :: Pointer image,InternalImageType :: Pointer maskingImage, IndexType pixels, IndexType seedPos, char * filename)
        {

            // Image size and spacing parameters
            unsigned long xExtent = size[0];
            unsigned long yExtent = size[1];
            unsigned long zExtent = size[2];
            unsigned long sourceImageSize[]  = { xExtent, yExtent, zExtent };


            double sourceImageSpacing[] = {image->GetSpacing()[0],image->GetSpacing()[1],image->GetSpacing()[2]};  // (pixels)
            //std::cout<< "sourceImageSpacing: " << sourceImageSpacing[0] << " x " << sourceImageSpacing[1] << " x " << sourceImageSpacing[2] << " mm^3" << std::endl;

            double sourceImageOrigin[] = {image->GetOrigin()[0],image->GetOrigin()[1],image->GetOrigin()[2]};  // (mm/pixels)
            //std::cout<< "sourceImageOrigin: " << sourceImageOrigin[0] << " x " << sourceImageOrigin[1] << " x " << sourceImageOrigin[2] << " mm^3" << std::endl;

            // Creates the sourceImage (but doesn't set the size or allocate memory)
            typename InternalImageType::Pointer sourceImage = InternalImageType::New();
            sourceImage->SetOrigin(sourceImageOrigin);
            sourceImage->SetSpacing(sourceImageSpacing);
            //std::cout<< "sourceImageSpacing: " << sourceImageSpacing[0] << " x " << sourceImageSpacing[1] << " x " << sourceImageSpacing[2] << " mm^3" << std::endl;

            //-----The following block allocates the sourceImage-----

            typename InternalImageType::SizeType sourceImageSizeObject;   // Create a size object native to the sourceImage type
            typename InternalImageType::RegionType largestPossibleRegion; // Create a region object native to the sourceImage type

            sourceImageSizeObject.SetSize( sourceImageSize );             // Set the size object to the array defined earlier
            largestPossibleRegion.SetSize( sourceImageSizeObject );       // Resize the region
            sourceImage->SetLargestPossibleRegion(largestPossibleRegion); // Set the largest legal region size (i.e. the size of the whole sourceImage) to what we just defined
            sourceImage->SetBufferedRegion( largestPossibleRegion );      // Set the buffered region
            sourceImage->SetRequestedRegion( largestPossibleRegion );     // Set the requested region
            sourceImage->Allocate();                                      // Now allocate memory for the sourceImage


            //std::cout << "New physical sourceImage allocated\n";


            //################################ Create ellipsoid in sourceImage ##################################

            // Ellipsoid spatial function typedef
            typedef itk::EllipsoidInteriorExteriorSpatialFunction<Dimension> EllipsoidFunctionType;
            typename EllipsoidFunctionType::Pointer spatialFunc = EllipsoidFunctionType::New();// Create an ellipsoid spatial function for the source image

            // Define and set the axes lengths for the ellipsoid
            typename EllipsoidFunctionType::InputType axes;
            axes[0] =  pixels[0];
            axes[1] =  pixels[1];
            axes[2] =  pixels[2];
            spatialFunc->SetAxes(axes);

            //ReaderWriter/ Define and set the center of the ellipsoid in physical space
            typename EllipsoidFunctionType::InputType center;
            center[0] = seedPos[0]* sourceImageSpacing[0];  // world coordinates(mm) -> pixels * spacing
            center[1] =seedPos[1]* sourceImageSpacing[1];  // world coordinates(mm) -> pixels * spacing
            center[2] = seedPos[2]* sourceImageSpacing[2];  // world coordinates(mm) -> pixels * spacing
            spatialFunc->SetCenter(center);

            double data[] = {0, 1, 0, 1, 0, 0, 0, 0, 1};
            vnl_matrix<double> orientations (data, 3, 3);

            // Set the orientations of the ellipsoids
            spatialFunc->SetOrientations(orientations);

            typedef  typename IndexType::IndexValueType   IndexValueType;

            IndexType seed;
            seed[0] = static_cast<IndexValueType>( seedPos[0] ); // (pixels)
            seed[1] = static_cast<IndexValueType>( seedPos[1] ); // (pixels)
            seed[2] = static_cast<IndexValueType>( seedPos[2]); // (pixels)

            itk::FloodFilledSpatialFunctionConditionalIterator<InternalImageType, EllipsoidFunctionType>
                    sfi = itk::FloodFilledSpatialFunctionConditionalIterator<InternalImageType,
                    EllipsoidFunctionType>(sourceImage, spatialFunc, seed);

            std::cout << "Seed position : " << seed << std::endl;

            // Iterate through the entire image and set interior pixels to 255
            int numInteriorPixels1 = 0;
            unsigned char interiorPixelValue = 255;
            // double intensity = 0;
            for(sfi.GoToBegin(); !sfi.IsAtEnd(); ++sfi) // sfi.GoToBegin()
            {
                //std::cout << "Index evaluated: " << sfi.GetIndex() << std::endl;
                sfi.Set(interiorPixelValue);
                ++numInteriorPixels1;
                //intensity += maskingImage->GetPixel(sfi.GetIndex());
            }
            //InternalImageType :: Pointer heighImage;
            //heighImage = HeightFunction(maskingImage);
            // ######################################### Mask image filter ##########################################

            typedef itk::MaskImageFilter< InternalImageType, InternalImageType > MaskFilterType;
            typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();

            maskFilter->SetInput(maskingImage);
            maskFilter->SetMaskImage(sourceImage);
            std::stringstream sss;

            sss << seed[0]<<"_"<<seed[1]<<"_"<<seed[2];
            std:: string strn = sss.str();
            string hei = filename + strn +"_mask.dcm";
            typedef itk::CastImageFilter< InternalImageType, OutputImageType > CastingFilterType1;
            CastingFilterType1::Pointer Castingfilter1 = CastingFilterType1::New();
            Castingfilter1->SetInput(maskFilter->GetOutput());
            try
            {
                toolsITK::writeDCMITKImage(Castingfilter1->GetOutput(),hei);

            }
            catch (itk::ExceptionObject & e)
            {
                std::cerr << "Exception caught !" << std::endl;
                std::cerr << e << std::endl;
                return 0;
            }

            return maskFilter->GetOutput();
        }

    public : int watershedSegmentation( InternalImageType :: Pointer image, float threshold,float level, IndexType index, char *outputfilename, char *filename)
        {
            //Watershed segmentation
            typedef itk::WatershedImageFilter<InternalImageType> WatershedFilterType;
            typename WatershedFilterType::Pointer watershedFilter=WatershedFilterType::New();

            //Gradient Magnitude Filter Output as Watershed Input
            watershedFilter->SetInput(image);
            //Threshold
            watershedFilter->SetThreshold(threshold);
            //Level
            watershedFilter->SetLevel(level);
            std::cout<<"Performing Watershed Segmentation..."<<std::endl;
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

            OutputPixelType value = watershedFilter->GetOutput()->GetPixel(index);
            std::cout<<"the tresholded intensity is "<< value <<endl;
            float newlevel=watershedFilter->GetLevel();;
            if (value == 1)
            {  while(value <= 1 && newlevel > 0 )
                {
                    newlevel = watershedFilter->GetLevel();
                    newlevel = newlevel - 0.01;
                    std::cout<<"new level value"<<newlevel<<endl;
                    std::cout<<"Performing Watershed Segmentation..."<<std::endl;
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
                    value = watershedFilter->GetOutput()->GetPixel(index);
                    std::cout<<"the tresholded intensity is "<< value <<endl;

                }
                std::cout<<"the tresholded intensity is "<< value <<endl;
            }
            std::cout<<"...Watershed segmentation Success(" << t.GetMean() << " s./" << t.GetMean()/60.0  << " min.)" <<std::endl;

            std::stringstream ss;
            ss << newlevel<<value;
            std:: string temp = ss.str();
            string waterShedfile = filename + temp + "watershed.nrrd";
            //WriteWatershedIntermediateimage(watershedFilter->GetOutput(),waterShedfile);
            typedef itk::CastImageFilter< WatershedImageType, OutputImageType > CastingFilterType2;
            CastingFilterType2::Pointer Castingfilter2 = CastingFilterType2::New();
            Castingfilter2->SetInput(watershedFilter->GetOutput());
            Castingfilter2->Update();

            //Seed Label
           // WriteDiacomImage(watershedFilter->GetOutput(),waterShedfile);

            try
            {
               toolsITK::writeDCMITKImage(Castingfilter2->GetOutput(),waterShedfile);
            }
            catch (itk::ExceptionObject & e)
            {
                std::cerr << "exception in file writer " << std::endl;
                std::cerr << e << std::endl;
                return 0;
            }
            
            OutputPixelType valueInten = watershedFilter->GetOutput()->GetPixel(index);

            //Label Threshold
            typedef itk::BinaryThresholdImageFilter<WatershedImageType, OutputImageType> ThresholdType;
            typename ThresholdType::Pointer thresholdFilter = ThresholdType::New();

            //Watershed Segmentation Output as Threshold Output
            thresholdFilter->SetInput(watershedFilter->GetOutput());
            thresholdFilter->SetInPlace(true);
            //Set Threshold filter limits
            thresholdFilter->SetLowerThreshold(valueInten);
            thresholdFilter->SetUpperThreshold(valueInten);
            //Set Output Values
            thresholdFilter->SetOutsideValue(0);
            thresholdFilter->SetInsideValue(1);
            //	std::cout<<"Extracting Seed Label..."<<std::endl;
            t.Start();
            //Apply Threshold
            thresholdFilter->Update();
            t.Stop();
            // Compute some statistics
            double sourceImageSpacing[] = {image->GetSpacing()[0],image->GetSpacing()[1],image->GetSpacing()[2]};  // (pixels)
            //std::cout<< "Image spacing: " << sourceImageSpacing[0] << " x " << sourceImageSpacing[1] << " x " << sourceImageSpacing[2] << " mm^3" << std::endl;

            getStatatistics(thresholdFilter->GetOutput(),sourceImageSpacing);

            try
            {
                toolsITK::writeDCMITKImage(thresholdFilter->GetOutput(),outputfilename);
            }
            catch (itk::ExceptionObject & e)
            {
                std::cerr << "exception in file writer " << std::endl;
                std::cerr << e << std::endl;
                return EXIT_FAILURE;
            }
            return EXIT_SUCCESS;
        }
    public : int** getCoordinates(IndexType seed, IndexType axex, double sourceImageSpacing[])
        {
            int** p = 0;
            p = new int*[8];

            for (int h = 0; h < 8; h++)
            {
                p[h] = new int[3];
            }
          //  std::cout<< "sourceImageSpacing: " << sourceImageSpacing[0] << " x " << sourceImageSpacing[1] << " x " << sourceImageSpacing[2] << " mm^3" << std::endl;

           axex[0] = axex[0]/(2*sourceImageSpacing[0]);          //pixels as axes is in mm
           axex[1] = axex[1]/(2*sourceImageSpacing[1]);
           axex[2] = axex[2]/(2*sourceImageSpacing[2]);

           //cout<< axex <<endl;
           // cout<< seed <<endl;
            p[0][0] = seed[0]-axex[0];                          //pixels as axes is in mm
            p[0][1] = seed[1]-axex[1];
            p[0][2] = seed[2]-axex[2];

            p[1][0] = seed[0]+axex[0];
            p[1][1] = seed[1]-axex[1];
            p[1][2] = seed[2]-axex[2];

            p[2][0] = seed[0]+axex[0];
            p[2][1] = seed[1]+axex[1];
            p[2][2] = seed[2]-axex[2];

            p[3][0] = seed[0]-axex[0];
            p[3][1] = seed[1]+axex[1];
            p[3][2] = seed[2]-axex[2];

            p[4][0] = seed[0]-axex[0];
            p[4][1] = seed[1]-axex[1];
            p[4][2] = seed[2]+axex[2];

            p[5][0] = seed[0]+axex[0];
            p[5][1] = seed[1]-axex[1];
            p[5][2] = seed[2]+axex[2];

            p[6][0] = seed[0]+axex[0];
            p[6][1] = seed[1]+axex[1];
            p[6][2] = seed[2]+axex[2];

            p[7][0] = seed[0]-axex[0];
            p[7][1] = seed[1]+axex[1];
            p[7][2] = seed[2]+axex[2];

            // for (int i = 0; i<8 ; i++)
            // {
            //     for (int j=0 ; j <3 ; j++)
            //     {
            //         cout<< p[i][j]<< " ";
            //     }
            //     cout<<endl;
            // }
            return p;

        }
    public: float** contrastRatio(InternalImageType :: Pointer image,InternalImageType :: Pointer maskingImage, IndexType pixels, int ** cubeCordinates, int length,char * filename)
        {

            // Image size and spacing parameters
            unsigned long xExtent = size[0];
            unsigned long yExtent = size[1];
            unsigned long zExtent = size[2];
            // unsigned long sourceImageSize[]  = { xExtent, yExtent, zExtent };
            // Now allocate memory for the sourceImage
            //std::cout << "New physical sourceImage allocated\n";
            //################################ Create ellipsoid in sourceImage ##################################

            typedef itk::MaskImageFilter< InternalImageType, InternalImageType > MaskFilterType;
            typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();

            typedef  typename IndexType::IndexValueType   IndexValueType;
            IndexType seedPos;

            float** vec = 0;
            vec = new float*[8];
            //std::cout<<sizeof(cubeCordinates)<<endl;
            for (int h = 0; h < length; h++)
            {
                vec[h] = new float[2];


            }

            unsigned long sourceImageSize[]  = { xExtent, yExtent, zExtent };

            double sourceImageSpacing[] = {image->GetSpacing()[0],image->GetSpacing()[1],image->GetSpacing()[2]};  // (pixels)
            //std::cout<< "sourceImageSpacing: " << sourceImageSpacing[0] << " x " << sourceImageSpacing[1] << " x " << sourceImageSpacing[2] << " mm^3" << std::endl;

            double sourceImageOrigin[] = {image->GetOrigin()[0],image->GetOrigin()[1],image->GetOrigin()[2]};  // (mm/pixels)
            //std::cout<< "sourceImageOrigin: " << sourceImageOrigin[0] << " x " << sourceImageOrigin[1] << " x " << sourceImageOrigin[2] << " mm^3" << std::endl;

            // Creates the sourceImage (but doesn't set the size or allocate memory)
            typename InternalImageType::Pointer sourceImage = InternalImageType::New();
            sourceImage->SetOrigin(sourceImageOrigin);
            sourceImage->SetSpacing(sourceImageSpacing);

            //std::cout<< "sourceImageSpacing: " << sourceImageSpacing[0] << " x " << sourceImageSpacing[1] << " x " << sourceImageSpacing[2] << " mm^3" << std::endl;

            //-----The following block allocates the sourceImage-----

            typename InternalImageType::SizeType sourceImageSizeObject;   // Create a size object native to the sourceImage type
            typename InternalImageType::RegionType largestPossibleRegion; // Create a region object native to the sourceImage type

            sourceImageSizeObject.SetSize( sourceImageSize );             // Set the size object to the array defined earlier
            largestPossibleRegion.SetSize( sourceImageSizeObject );       // Resize the region
            sourceImage->SetLargestPossibleRegion(largestPossibleRegion); // Set the largest legal region size (i.e. the size of the whole sourceImage) to what we just defined
            sourceImage->SetBufferedRegion( largestPossibleRegion );      // Set the buffered region
            sourceImage->SetRequestedRegion( largestPossibleRegion );     // Set the requested region
            sourceImage->Allocate(); // Now allocate memory for the sourceImage
            // getStatatistics(sourceImage);
            // Ellipsoid spatial function typedef
            typedef itk::EllipsoidInteriorExteriorSpatialFunction<Dimension> EllipsoidFunctionType;
            typename EllipsoidFunctionType::Pointer spatialFunc = EllipsoidFunctionType::New();// Create an ellipsoid spatial function for the source image

            // Define and set the axes lengths for the ellipsoid
            typename EllipsoidFunctionType::InputType axes;
            axes[0] =  pixels[0]//*sourceImageSpacing[0] ;
            axes[1] =  pixels[1]//*sourceImageSpacing[1] ;
            axes[2] =  pixels[2]//*sourceImageSpacing[2] ;
            spatialFunc->SetAxes(axes);

           
            for (int i = 0; i< length ; i++) {           

            
                //ReaderWriter/ Define and set the center of the ellipsoid in physical space
                typename EllipsoidFunctionType::InputType center;
                double sigm = 0, mean=0;
                seedPos[0]= static_cast<IndexValueType>(cubeCordinates[i][0]);
                seedPos[1]= static_cast<IndexValueType>(cubeCordinates[i][1]);
                seedPos[2]= static_cast<IndexValueType>(cubeCordinates[i][2]);
                if(seedPos[0]== xExtent)
                    seedPos[0] = seedPos[0]-pixels[0]-1;
                if(seedPos[2]< 0)
                    seedPos[2] = pixels[0];
                center[0] = seedPos[0]* sourceImageSpacing[0];  // world coordinates(mm) -> pixels * spacing
                center[1] = seedPos[1]* sourceImageSpacing[1];  // world coordinates(mm) -> pixels * spacing
                center[2] = seedPos[2]* sourceImageSpacing[2];  // world coordinates(mm) -> pixels * spacing
                spatialFunc->SetCenter(center);

                double data[] = {0, 1, 0, 1, 0, 0, 0, 0, 1};
                vnl_matrix<double> orientations (data, 3, 3);

                // Set the orientations of the ellipsoids
                spatialFunc->SetOrientations(orientations);

                itk::FloodFilledSpatialFunctionConditionalIterator<InternalImageType, EllipsoidFunctionType>
                        sfi = itk::FloodFilledSpatialFunctionConditionalIterator<InternalImageType,
                        EllipsoidFunctionType>(sourceImage, spatialFunc, seedPos);

                //  std::cout << "Seed position : " << seedPos << std::endl;

                // Iterate through the entire image and set interior pixels to 255
                int numInteriorPixels1 = 0;
                unsigned char interiorPixelValue = 255;
                double intensity = 0;
                for(sfi.GoToBegin(); !sfi.IsAtEnd(); ++sfi) // sfi.GoToBegin()
                {
                    //std::cout << "Index evaluated: " << sfi.GetIndex() << std::endl;
                    sfi.Set(interiorPixelValue);
                    ++numInteriorPixels1;
                    intensity += maskingImage->GetPixel(sfi.GetIndex());
                }
                //std::cout << "No. voxels included in sphere mean : " << numInteriorPixels1 << std::endl;
                // std::cout << "Sum of intensity : " << intensity << std::endl;
                mean = intensity/ numInteriorPixels1;
                std::cout << "Mean of the region " << i << " is : "<< mean << std::endl;
                double variance =0;
                numInteriorPixels1 = 0;
                for(sfi.GoToBegin(); !sfi.IsAtEnd(); ++sfi) // sfi.GoToBegin()
                {
                    //                //std::cout << "Index evaluated: " << sfi.GetIndex() << std::endl;
                    //                sfi.Set(interiorPixelValue);
                    ++numInteriorPixels1;
                    variance += (maskingImage->GetPixel(sfi.GetIndex())-mean)*(maskingImage->GetPixel(sfi.GetIndex())-mean);
                }

                sigm = sqrt(variance / numInteriorPixels1);
                std::cout << "Standard deviation of the region " << i << " is:  "<<sigm << std::endl;
               // std::cout << "No. voxels included in sphere variance : " << numInteriorPixels1 << std::endl;

                vec[i][0] = mean;
                vec[i][1] = sigm;

                maskFilter->SetInput(maskingImage);
                maskFilter->SetMaskImage(sourceImage);
                           

        }
            std::stringstream sss;
            sss << "co-ordinate" ;//<< seedPos[0]<<"_"<<seedPos[1]<<"_"<<seedPos[2];
            std:: string strn = sss.str();
            string hei = filename + strn +"_mask.dcm";

            typedef itk::CastImageFilter< InternalImageType, OutputImageType > CastingFilterType1;
            CastingFilterType1::Pointer Castingfilter1 = CastingFilterType1::New();
            Castingfilter1->SetInput(maskFilter->GetOutput());
            try
            {
                toolsITK::writeDCMITKImage(Castingfilter1->GetOutput(),hei);

            }
            catch (itk::ExceptionObject & e)
            {
                std::cerr << "Exception caught !" << std::endl;
                std::cerr << e << std::endl;
                return 0;
            }
            return vec;

        }
    public: float** getGTMask (InternalImageType :: Pointer image, InternalImageType :: Pointer maskingImage, char * filename ){

           // std::cout << "Here";
            typedef itk::MaskImageFilter< InternalImageType, InternalImageType > MaskFilterType;
            typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();
           // std::cout << "Here";
            maskFilter->SetInput(image);
            maskFilter->SetMaskImage(maskingImage);
            std::stringstream sss;

            sss << "_";
            std:: string strn = sss.str();

            string hei = filename +strn + "maskGT.dcm";

            typedef itk::CastImageFilter< InternalImageType, OutputImageType > CastingFilterType1;
            CastingFilterType1::Pointer Castingfilter1 = CastingFilterType1::New();
            Castingfilter1->SetInput(maskFilter->GetOutput());
            try
            {
                toolsITK::writeDCMITKImage(Castingfilter1->GetOutput(),hei);

            }
            catch (itk::ExceptionObject & e)
            {
                std::cerr << "Exception caught !" << std::endl;
                std::cerr << e << std::endl;
                return 0;
            }
            
        }

    public: float** getCnrSNR(maskingImage)
    {
            float intensity = 0.0;
            int numInteriorPixels1 = 0;
            //InternalImageType::Pointer im = maskingImage;
            itk::ImageRegionIterator<InternalImageType>  it( maskingImage, maskingImage->GetLargestPossibleRegion() );
            for (it.GoToBegin(); !it.IsAtEnd(); ++it)
            {
                if(maskingImage->GetPixel(it.GetIndex()))
                {
                    numInteriorPixels1++;
                    intensity += maskFilter->GetOutput()->GetPixel(it.GetIndex());
                }
            }
            float mean = 0.0, sigma = 0.0;
            mean = intensity/ numInteriorPixels1;
            std::cout << "Mean of the region is "<<mean <<" "<< std::endl;
            double variance =0;
            // numInteriorPixels1 = 0;
            for(it.GoToBegin(); !it.IsAtEnd(); ++it) // sfi.GoToBegin()
            {
                if(maskingImage->GetPixel(it.GetIndex()))
                {
                    variance += ( maskFilter->GetOutput()->GetPixel(it.GetIndex())-mean)*(maskFilter->GetOutput()->GetPixel(it.GetIndex())-mean);
                }
            }
            std::cout << "variance of the region is "<<variance << std::endl;
            sigma= sqrt(variance / numInteriorPixels1);
            std::cout << "Sigma of the region is "<<sigma << std::endl;

           float** vec = 0;
           vec = new float*[1];
           //std::cout<<sizeof(cubeCordinates)<<endl;
           for (int h = 0; h < 1; h++)
           {
               vec[h] = new float[2];
           }
           vec[0][0] = mean;
           vec[0][1] = sigma;
            return vec;
}
    };


    int main(int argc, char *argv[]) {


        //cout<< argv[10]<<endl;
        if(argc!=11)
        {
            std::cout<<"Input_image Output_image SeedX SeedY SeedZ"<<std::endl;
            return 0;
        }
        

        itk::TimeProbe clock;
        clock.Start();
        Segmentation watershedseg;
        Segmentation groundTruth;
        
        int nei =5;                    //mm diameter for each background  pixel     
        float conductance = 3;  
        float level=0.60;               // Controls watershed depth. Percentage(0.0-1.0) of the maximum depth in the input image.
        float threshold = 0.005;        // 0.1 - controls the lower thresholding of the input. Percentage(0.0-1.0) of the maximum depth in the input image.
                                        // (A rule of thumb is to set the Threshold to be about 1 / 100 of the Level.)
        int iterations =5;              // Typical 5. more interactions will smooth further but will increase computing time.
      //float sigma = atof(argv[10]);
        float sigma =0.1; 

                // 1 - Around 3.0
        char * name = argv[9];


        IndexType seedPoint;
        seedPoint[0] = atoi( argv[3] ); // (pixels)
        seedPoint[1] = atoi( argv[4] ); // (pixels)
        seedPoint[2] = atoi( argv[5] ); // (pixels)

        // int** seed = 0;
        // seed = new int*[5];
        // for (int h = 0; h < 1; h++)
        // {
        //     seed[h] = new int[3];
        // }

        // seed[0][0] = seedPoint[0];   // (pixels)
        // seed[0][1] = seedPoint[1];   // (pixels)
        // seed[0][2] = seedPoint[2];   // (pixels)

        IndexType axesbackground;
        axesbackground[0] = atoi( argv[6] );      // (mm)
        axesbackground[1] = atoi( argv[7] );      // (mm)
        axesbackground[2] = atoi( argv[8] );      // (mm)

        IndexType axes;
        axes[0] = 30;      // (mm)
        axes[1] = 30;      // (mm)
        axes[2] = 30;      // (mm)


        
        //std::cout<< "cubeCordinates:: "<< cubeCordinates[0][1]<<"  " <<cubeCordinates[1][2]<<endl;

        // std::cout<< "cubeCordinates:: "<< point<<endl;
        IndexType neighbours;
        neighbours[0] = nei;            // (pixels)
        neighbours[1] = nei;            // (pixels)
        neighbours[2] = nei;            // (pixels)
        int count =0;
        ofstream myfile;
        
        //Perform algorithm

        watershedseg.image = watershedseg.ReadImage(argv[1]);           //load the image

        groundTruth.image = groundTruth.ReadImage(argv[10]);            //load the ground truth

        InternalImageType :: Pointer  maskImage = watershedseg.CreateMask(watershedseg.image,watershedseg.image,axesbackground,seedPoint, name);  //Truncating the image to get a mask od 50mm

        double sourceImageSpacing[] = {watershedseg.image->GetSpacing()[0],watershedseg.image->GetSpacing()[1],watershedseg.image->GetSpacing()[2]}; //get the spacing of the image

        int** cubeCordinates= watershedseg.getCoordinates(seedPoint, axes, sourceImageSpacing);     //generate the 8 co-ordinates around the mask for background information
        //std::cout<< "cubeCordinates:: "<< cubeCordinates[0][0]<<"  " <<cubeCordinates[0][1]<<"  " <<cubeCordinates[0][2]<<endl;

       float ** sigMuStd = watershedseg.contrastRatioGT(watershedseg.image,groundTruth.image,name);
        sss << "_";
        std:: string strn = sss.str();
        string hei = name +strn + "maskGT.dcm";
        groundTruth.image = groundTruth.ReadImage(hei); 


            
        float ** bgMuSTD = watershedseg.contrastRatio(watershedseg.image,watershedseg.image,neighbours,cubeCordinates,8,name);

        float bgMean= 0, bgSigma= 0;

        for(int i =0;i<8;i++){
            bgMean += bgMuSTD[i][0];
            bgSigma += bgMuSTD[i][1];
        }

        bgMean = bgMean/8;
        bgSigma = bgSigma/8;


        float cnr = abs(sigMuStd[0][0] - bgMean)/sqrt(sigMuStd[0][1]*sigMuStd[0][1] + bgSigma*bgSigma);
        float snr = sigMuStd[0][0]/bgSigma; 
        myfile.open ("ratio.csv",std::ios_base::app);
        myfile<<"\n" ;
        myfile << count<<","<<cnr<<"\n";
        myfile.close();

        std::cout << "###############################---------"<<sigMuStd[0][0]<<" , " << sigMuStd[0][1] <<", "  <<bgMean<<" , "<< bgSigma<<std::endl;
        std::cout << "###############################----CNR-----"<<cnr<<endl;
        std::cout << "###############################----SNR-----"<<snr<<endl;

        float cnrNew = 0;
        float ** muSigNew = 0;
        float snrNew = 0;

        InternalImageType :: Pointer smoothimage =   watershedseg.anisotropicDiffusion(maskImage,iterations,conductance);
         InternalImageType :: Pointer gtSmoothImage =   watershedseg.anisotropicDiffusion(groundTruth.image,iterations,conductance);
        
        InternalImageType :: Pointer  gaussianSmoothimage = watershedseg.HeightFunction(smoothimage,sigma);

        //do {

       
           float bgMeanNew= 0, bgSigmaNew= 0;
           float ** bgMuSTDNew = watershedseg.contrastRatio(watershedseg.image,gaussianSmoothimage,neighbours,cubeCordinates,8,name);
            for(int i =0;i<8;i++)
            {
                bgMeanNew += bgMuSTDNew[i][0];
                bgSigmaNew += bgMuSTDNew[i][1];
            }
            bgMeanNew = bgMeanNew/8;
            avSigmaNew = avSigmaNew/8;
            //muSigNew = watershedseg.contrastRatio(watershedseg.image,smoothimage,neighbours,seed,1,name);
            cnrNew = abs(sigMuStd[0][0] - avMeanNew)/sqrt(sigMuStd[0][1]*sigMuStd[0][1] + avSigmaNew*avSigmaNew);
            snrNew   = sigMuStd[0][0]/avSigmaNew;
            std::cout<<"Iteration ::"<< count<<endl;
            std::cout << "###############################---------"<<muSigNew[0][0]<<" , " << muSigNew[0][1] <<", "  <<avMeanNew<<" , "<< avSigmaNew<<std::endl;
            std::cout << "###############################----CNR-----"<<cnrNew<<endl;

            std::cout << "###############################----SNR-----"<<snrNew<<endl;

        //    if((cnrNew-cnr)< 0.1 )
        //    {    cnr = cnrNew;
        //        break;
        //    }
        //    cnr = cnrNew;
        //     smoothimage =   watershedseg.anisotropicDiffusion(smoothimage,iterations,conductance);;
        //     count += 5;
        // }while(count <= 35);
        
        // myfile.open ("ratio_iterations.csv",std::ios_base::app);
        // myfile << count<<","<<cnrNew<<"\n";
        // myfile.close();



        
        
        // do {


        //     float avMeanNew= 0, avSigmaNew= 0;
        //     float ** bgMuSTDNew = watershedseg.contrastRatio(watershedseg.image,gaussianSmoothimage,neighbours,cubeCordinates,8,name);
        //     for(int i =0;i<8;i++)
        //     {
        //         avMeanNew += bgMuSTDNew[i][0];
        //         avSigmaNew += bgMuSTDNew[i][1];
        //     }
        //     avMeanNew = avMeanNew/8;
        //     avSigmaNew = avSigmaNew/8;
        //     muSigNew = watershedseg.contrastRatio(watershedseg.image,gaussianSmoothimage,neighbours,seed,1,name);
        //     cnrNew = abs(muSigNew[0][0] - avMeanNew)/sqrt(muSigNew[0][1]*muSigNew[0][1] + avSigmaNew*avSigmaNew);
        //     snrNew   = muSigNew[0][0]/muSigNew[0][1];
        //     std::cout<<"Sigma ::"<< sigma<<endl;
        //     std::cout << "###############################---------"<<muSigNew[0][0]<<" , " << muSigNew[0][1] <<", "  <<avMeanNew<<" , "<< avSigmaNew<<std::endl;
        //     std::cout << "###############################----CNR-----"<<cnrNew<<endl;

        //     std::cout << "###############################----SNR-----"<<snrNew<<endl;

        //     if(cnrNew < cnr )
        //     {
        //         gaussianSmoothimage = watershedseg.Smoothing(smoothimage,sigma);
        //         break;
        //     }
        //     cnr = cnrNew;

        //     sigma += 0.3;

        // }while(sigma <= 1.5);

        myfile.open ("ratio_sigma.csv",std::ios_base::app);
        myfile << sigma<<","<<cnr<<"\n";
        myfile.close();

        
        watershedseg.watershedSegmentation(maskImage, threshold, level, seedPoint, argv[2],name);

        // print computation time at screen
        clock.Stop();
        std::cout << "###############################" << std::endl;
        std::cout << "SegmentationABUS" << std::endl;
        std::cout << "Total time: " << clock.GetTotal()/60.0  << " min." << std::endl;
        std::cout << "Total time: " << clock.GetTotal()       << " s." << std::endl;
        std::cout << "###############################" << std::endl;

        return 0;
    }



