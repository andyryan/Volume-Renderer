#include "colour.h"
#include "Image.h"
#include "Volume.h"

#include <stdio.h>
#include <string>

int main( int argc, const char* argv[] )
{
	/*Image image;
	std::string filename = "lena256.bmp";
	std::ofstream imageFile;
	imageFile.open("output.bmp", std::ios_base::binary);
	image.readBitmap(filename);
	image.writeBitmap(imageFile);*/
	Volume volume;
	std::string fileName = "foot.raw";
	volume.readRaw(fileName, 256, 256, 256, 1);
	std::ofstream volumeFile;
	volumeFile.open("output.raw", std::ios_base::binary);


}