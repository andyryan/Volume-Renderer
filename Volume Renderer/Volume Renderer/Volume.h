#pragma once

#include <cmath>
#include <string>
#include <fstream>
#include "Colour.h"
#include "Image.h"
#include "Ray.h"

struct HitRecord{
	float t;
	Vector3 normal;
	Colour colour;
};

class Volume
{
public:
	Volume(void);
	~Volume(void);
	Volume(int width, int height, int depth, int channels);
	void clearBuffer();
	void allocateBuffer();
	void setBoundingBox();

	void writePVM(std::ostream& out);
	void readPVM(std::string file_name);
	void writeRaw(std::ostream& out);
	void readRaw(std::string file_name, int width, int height, int planes, int channels);

	Image* maximumIntensityProjection(int width, int height, float theta, float phi = 0);
	Volume & backProject(Image & sinogram);

	bool hit(const Ray& r, float tmin, float tmax, float time, HitRecord& hit) const;


private:
	Colour ***buffer;
	int boundingBox[6];
	int _width, _height, _planes, _channels;
};

