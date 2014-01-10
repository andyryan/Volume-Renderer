#include "Image.h"
#include <iostream>
#include <Windows.h>

Image::Image(){
	raster = NULL;
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

Image::Image(std::string fileName){
	readPPM(fileName);
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
	unsigned int ired, igreen, iblue;
	unsigned char cred, cgreen, cblue;

	for (j = 0; j< _height; j++){
		for (i = 0; i < _width; i++){
			if (isRGB() || isRGBA()){
				ired = (unsigned int)(255*raster[i][j][0]);
				igreen = (unsigned int)(255*raster[i][j][1]);
				iblue = (unsigned int)(255*raster[i][j][2]);
			}
			if (isGreyscale()){
				ired = (unsigned int)(255*raster[i][j][0]);
				igreen = (unsigned int)(255*raster[i][j][0]);
				iblue = (unsigned int)(255*raster[i][j][0]);
			}
			if (ired > 255) ired = 255;
			if (igreen > 255) igreen = 255;
			if (iblue > 255) iblue = 255;
			cred = (unsigned char)(ired);
			cgreen = (unsigned char)(igreen);
			cblue = (unsigned char)(iblue);
			out.put(cred);
			out.put(cgreen);
			out.put(cblue);
		}
	}
}

void Image::readPPM(std::string fileName){
//open stream to file
	std::ifstream in;
	in.open (fileName.c_str(), std::ios_base::binary);
	if(!in.is_open()){
		std::cerr << " ERROR -- Couldn't open file \'" << fileName << "\'.\n";
		exit(-1);
	}

	char ch, type;
	char red, green, blue;
	int i, j, cols, rows;
	int num;

	//get header info
	in.get(ch);
	in.get(type);
	in >> cols >> rows >> num;

	_width = cols;
	_height = rows;
	_channels = 3;

	//allocate raster
	if (raster != NULL){
		for (i = 0; i < _width; i++)
		if (raster[i] != NULL)
			delete raster[i];
		delete raster;
	}
	raster = new Colour*[_width];
	for (i = 0; i < _width; i++)
		raster[i] = new Colour[_height];

	//cleanup newline
	in.get(ch);

	//store ppm values
	for ( j = 0; j < _height; j++)
		for ( i = 0; i < _width; i++){
			in.get(red);
			in.get(green);
			in.get(blue);

			float values[3] = {(float)((unsigned char)red)/255.0,
							   (float)((unsigned char)green)/255.0,
							   (float)((unsigned char)blue)/255.0};

			raster[i][j] = Colour(3, values);
		}
}

void Image::readBitmap(std::string fileName){
	FILE *file;
	BYTE *raw;
	BITMAPFILEHEADER bfh;
	BITMAPINFOHEADER bih;
	file = fopen(fileName.c_str(), "rb");
	if (file == NULL)
		return;
	fread(&bfh, sizeof(BITMAPFILEHEADER), 1,file);
	if (bfh.bfType != 0x4D42)
	{
		fclose(file);
		return;
	}

	fread(&bih, sizeof(BITMAPINFOHEADER),1,file);

	_width = bih.biWidth;
	_height = bih.biHeight;
	_channels = bih.biBitCount/8;
	//isFourier = false;
	//normalized = false;
	int i;
	
	if (raster != NULL){
		for (i = 0; i < _width; i++)
		if (raster[i] != NULL)
			delete raster[i];
		delete raster;
	}
	raster = new Colour*[_width];
	for (i = 0; i < _width; i++)
		raster[i] = new Colour[_height];

	int bufferSize = _channels*_width*_height;
	//deletePixels();
	//pixel = (double*)malloc(buffersize);
	raw = new BYTE[bufferSize];
	fseek(file, bfh.bfOffBits,SEEK_SET);
	fread(raw, 1, bufferSize, file);
	float *values = new float[_channels];
	for (int j = _height-1; j >=0; j--)
	{
		for (int i = 0; i < _width; i++)
		{
			int index = (i*_channels)+(j*_width*_channels);
			for (int k = 0; k < _channels; k++)
				values[(_channels-1) - k] = (float)raw[index+k]/255.0f;
			Colour c(_channels, values);
			set(i, (_height-1)-j, c);
		}
	}

	fclose(file);
	//normalize();
	//std::cout << "Image "<< fileName << " input correctly \n";
	//std::cout << "Width: "<< dim[1] << "\n";
	//std::cout << "Height: " << dim[2] << "\n";
	delete[] raw;
	delete values;
}

void Image::writeBitmap(std::ostream& out){
	if (isCOMP()) throw "Incompatible image for .bmp";
	//if (_channels == 1)
		//toColour();
	BYTE *raw;

	BITMAPFILEHEADER bfh;
	BITMAPINFOHEADER bih;
	bfh.bfType = 0x4D42;
	bfh.bfSize = _channels*_width*_height+sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER);
	bfh.bfReserved1 = 0;
	bfh.bfReserved2 = 0;
	bfh.bfOffBits= sizeof(BITMAPINFOHEADER)+sizeof(BITMAPFILEHEADER);

	bih.biSize = sizeof(BITMAPINFOHEADER);
	bih.biWidth = _width;
	bih.biHeight = _height;
	bih.biPlanes = 1;
	bih.biBitCount = 24;
	bih.biSizeImage = _width*_height*_channels;
	bih.biXPelsPerMeter = 2400;
	bih.biYPelsPerMeter = 2400;
	if (_channels > 1){
		bih.biClrUsed = 0;
		bih.biClrImportant =0;
	}
	else {
		bih.biClrUsed = 256;
		bih.biClrImportant =256;
	}
	bih.biCompression = BI_RGB;

	out.write((char*)&bfh, sizeof(bfh));
    out.write((char*)&bih, sizeof(bih));
	out.seekp(bfh.bfOffBits);
	//fwrite(&bfh, sizeof(BITMAPFILEHEADER),1,out);
	//fwrite(&bih, sizeof(BITMAPINFOHEADER),1,output);
	//fseek(output, bfh.bfOffBits,SEEK_SET);
	
	int bufferSize = _width*_height*_channels;
	raw = new BYTE[bufferSize];
	for (int j = _height-1; j >=0; j--)
	{
		for (int i = 0; i < _width; i++)
		{
			int index = (i*_channels)+(j*_width*_channels);
			for (int k = 0; k < _channels; k++)
			{
				unsigned int c = ((unsigned int)(raster[i][(_height-1)-j][(_channels-1)-k]*255));
				if (c > 255) c = 255;
				raw[index+k] = (unsigned char)c;
			}
		}
	}
	out.write((char*)raw, bufferSize);
	

	//std::cout << "Image "<< fileName << " output correctly \n";
	//std::cout << "Width: "<< dim[1] << "\n";
	//std::cout << "Height: " << dim[2] << "\n";
	delete[] raw;
	//fclose(output);

}