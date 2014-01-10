#include "Volume.h"
#include <iostream>
#include <Windows.h>


Volume::Volume(void)
{
	buffer = NULL;
}


Volume::~Volume(void)
{
}

void Volume::readPVM(std::string fileName){
//open stream to file
	std::ifstream in;
	in.open (fileName.c_str(), std::ios_base::binary);
	if(!in.is_open()){
		std::cerr << " ERROR -- Couldn't open file \'" << fileName << "\'.\n";
		exit(-1);
	}

	char ch;
	char type[7];
	char red, green, blue;
	int i, j, cols, rows, planes;
	int num;

	//get header info
	in.get(type, 7);
	in.get(ch);
	in >> cols >> rows >> planes;

	_width = cols;
	_height = rows;
	_planes = planes;
	_channels = 3;

	//allocate raster
	if (buffer != NULL){
		for (i = 0; i < _width; i++){
			if (buffer[i] != NULL){
				for (j = 0; j < _height; i++)
					if (buffer[i][j] != NULL)
						delete buffer[i][j];
				delete buffer[i];
			}
		}
		delete buffer;
	}
	buffer = new Colour**[_width];
	for (i = 0; i < _width; i++){
		buffer[i] = new Colour*[_height];
		for (j = 0; j < _height; j++)
			buffer[i][j] = new Colour[_planes];
	}

	//cleanup newline
	in.get(ch);

	//store ppm values
	char * colourCh = new char[_channels];
	float * colourF = new float[_channels]; 
	for (int z = 0; z < _planes; z++)
		for ( int y = 0; y < _height; y++)
			for ( int x = 0; x < _width; x++){
				for (int c = 0; c < _channels; c++){
					in.get(colourCh[c]);
					colourF[c] = (float)((unsigned char)red)/255.0;
				}
				buffer[x][y][z] = Colour(_channels, colourF);
			}
	delete colourCh;
	delete colourF;
}

void Volume::readRaw(std::string fileName, int width, int height, int planes, int channels){
//open stream to file
	std::ifstream in;
	in.open (fileName.c_str(), std::ios_base::binary);
	if(!in.is_open()){
		std::cerr << " ERROR -- Couldn't open file \'" << fileName << "\'.\n";
		exit(-1);
	}

	//allocate raster
	if (buffer != NULL){
		for (int i = 0; i < _width; i++){
			if (buffer[i] != NULL){
				for (int j = 0; j < _height; i++)
					if (buffer[i][j] != NULL)
						delete buffer[i][j];
				delete buffer[i];
			}
		}
		delete buffer;
	}

	_width = width;
	_height = height;
	_planes = planes;
	_channels = channels;

	buffer = new Colour**[_width];
	for (int i = 0; i < _width; i++){
		buffer[i] = new Colour*[_height];
		for (int j = 0; j < _height; j++)
			buffer[i][j] = new Colour[_planes];
	}
	//store ppm values
	char * colourCh = new char[_channels];
	float * colourF = new float[_channels]; 
	for (int z = 0; z < _planes; z++)
		for ( int y = 0; y < _height; y++)
			for ( int x = 0; x < _width; x++){
				for (int c = 0; c < _channels; c++){
					in.get(colourCh[c]);
					colourF[c] = (float)((unsigned char)colourCh[c])/255.0;
				}
				buffer[x][y][z] = Colour(_channels, colourF);
			}
	delete colourCh;
	delete colourF;
}

void Volume::writeRaw(std::ostream& out){
	for (int z = 0; z < _planes; z++)
		for ( int y = 0; y < _height; y++)
			for ( int x = 0; x < _width; x++){
				for (int c = 0; c < _channels; c++){
					out.put(buffer[x][y][z][c]);
				}
			}
}

Image* Volume::maximumIntensityProjection(int width, int height, float theta, float phi){
	Image *output= new Image(width, height, _channels);


	return output;
}
