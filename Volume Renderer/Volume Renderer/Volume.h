#pragma once

#include <cmath>
#include <string>
#include <fstream>
#include "Colour.h"
#include "Image.h"

class Volume
{
public:
	Volume(void);
	~Volume(void);

	void writePVM(std::ostream& out);
	void readPVM(std::string file_name);
	void writeRaw(std::ostream& out);
	void readRaw(std::string file_name, int width, int height, int planes, int channels);

	Image* maximumIntensityProjection(int width, int height, float theta, float phi = 0);
	Volume & backProject(Image & sinogram);


private:
	Colour ***buffer;
	int _width, _height, _planes, _channels;
};

