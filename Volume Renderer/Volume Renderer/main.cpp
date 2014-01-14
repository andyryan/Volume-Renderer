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

	HitRecord rec;
	bool is_a_hit;
	float tmax;
	Vector3 dir(0,0,1);

	Image output(256, 256, 1);
	for (int j = 0; j < 256; j++)
		for (int i = 0; i < 256; i++){
			tmax = 100000.0f;
			is_a_hit = false;
			Ray r(Vector3(i, j, 0), dir);
			if (volume.hit(r, 0.00001f, tmax, 0, rec)){
				tmax = rec.t;
				is_a_hit = true;
			}
			if (is_a_hit)
				output.set(i, j, rec.colour);
			else{
				float c[] = {0, 0, 0};
				output.set(i, j, Colour(1, c));
			}
		}

	std::ofstream imageFile;
	imageFile.open("projectionOutput.ppm", std::ios::trunc);
	output.writePPM(imageFile);
	imageFile.close();
	volumeFile.close();

}