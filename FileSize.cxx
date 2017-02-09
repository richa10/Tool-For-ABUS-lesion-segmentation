
#include "utility.h"

int main(int argc, char *argv[])
{
	//cout << argc << endl;
	if (argc != 5) {
		cout << "movingImage fixedImage patient file" << std::endl;
		return 0;
	}
	itk::TimeProbe clock;
	clock.Start();
	char * patient = argv[3];
	char * file = argv[4];

	InternalImageType::Pointer  movingImage, fixedImage;
	myfile2.open(file, ios::out | ios::app);
	myfile2 << patient << ",";
	ReadImage(argv[1]);
	ReadImage(argv[2]);
	myfile2 << endl;
}

