#include "colour.h"
#include "Image.h"

#include <stdio.h>
#include <string>

int main( int argc, const char* argv[] )
{
	Image image;
	std::string filename = "lena256.bmp";
	std::ofstream imageFile;
	imageFile.open("output.bmp", std::ios_base::binary);
	image.readBitmap(filename);
	image.writeBitmap(imageFile);

}