#include "Image.h"
#include <iostream>

Image::Image(){

}

Image::Image(int width, int height, int channels){
	_width = width;
	_height = height;
	_channels = channels;
	raster = new Colour*[_width];
	for (int i = 0; i < _width; i++){
		raster[i] = new Colour[_height];
	}
}

Image::Image(int width, int height, Colour background){
	_width = width;
	_height = height;
	_channels = background.nbChannels();
	raster = new Colour*[_width];
	for (int i = 0; i < _width; i++){
		raster[i] = new Colour[_height];
		for (int j = 0; j < _height; j++){
			raster[i][j] = background;
		}
	}
}

bool Image::set(int x, int y, const Colour& colour){
	if (x < 0 || x > _width) throw "x out of bounds";
	if (y < 0 || y > _width) throw "y out of bounds";
	if (colour.nbChannels() != _channels) throw "incompatible colour";
	raster[x][y] = colour;
	return true;
}

Colour Image::get(int x, int y){
	if (x < 0 || x > _width) throw "x out of bounds";
	if (y < 0 || y > _width) throw "y out of bounds";
	return raster[x][y];
}

void Image::writePPM(std::ostream& out){

	if (isCOMP()) throw "incompatible Image";
	//output header
	out << "P6\n";
	out << _width << ' ' << _height << '\n';
	out << "255\n";
	
	int i, j;
	unsigned int colourInt;
	unsigned char colourChar;

	for (i = _height-1; i>=0; i--){
		for (j = 0; j< _width; j++){
			for (int k = 0; k < _channels; k++){
				colourInt = (unsigned int)(255*raster[i][j][k]);
				if (colourInt > 255) colourInt = 255;
				colourChar = (unsigned char)colourInt;
				out.put(colourChar);
			}
		}
	}
}

void Image::readPPM(std::string fileName){
	std::ifstream in;
	in.open(fileName.c_str());
	if (!in.is_open()){
		std::cerr << "Error -- couldn't open file \'" << fileName << "\'. \n";
	}



}