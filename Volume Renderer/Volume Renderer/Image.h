#pragma once

#include <cmath>
#include <string>
#include <fstream>
#include "Colour.h"

class Image
{
public:
	Image();
	Image(int width, int height);
	Image(int width, int height, Colour background);
	Image(int width, int height, int nbChannels);

	bool set (int x, int y, const Colour& colour);
	Colour get(int x, int y);
	void writePPM(std::ostream& out);
	void readPPM(std::string file_name);
	void writeBitmap(std::ostream& out);
	void readBitmap(std::string file_name);
	void writeBin(std::ostream& out);
	void readBin(std::string file_name);

	bool isRGB();
	bool isRGBA();
	bool isGreyscale();
	bool isCOMPRGB();
	bool isCOMPGREY();

	void initWith1();
	void initWith0();

	Image& resize(int w, int h);
	Image& toGreyscale();
	Image& toColour();
	Image& toComplexColour();
	Image& toComplexGrey();


	Image& operator=(const Image & right_op);
	Image& operator+=(const Image& right_op);
	Image& operator*=(float s);
	Image& operator/=(float s);
	Image& operator*=(const Image& right_op);
	Image& operator/=(const Image& right_op);

	Image& fft(int direction);
	Image& fft(int direction, int axis);
	Image& fftShiftH();
	Image& fftShiftV();
	Image& fftShift();

	void ramLakFilter();
	void butterworthFilter();
	void sheppLoganFilter();
	void sheppLoganFilter3D();
	void cosFilter();
	void hannFilter();
	void hammingFilter();


	int width(){return _width;}
	int height(){return _height;}
	int channels();

	bool isValid();
	void clamp();
	Image& normalize();


private:
	Colour **raster;
	int _width, _height;
};

