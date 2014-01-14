#include "Volume.h"
#include <iostream>


Volume::Volume(void)
{
	buffer = NULL;
	boundingBox[0] = 1000000;
	boundingBox[1] = 0;
	boundingBox[2] = 1000000;
	boundingBox[3] = 0;
	boundingBox[4] = 1000000;
	boundingBox[5] = 0;
}

Volume::Volume(int width, int height, int depth, int channels){
	_width = width;
	_height = height;
	_planes = depth;
	buffer = new Colour**[_width];
	for (int i = 0; i < _width; i++){
		buffer[i] = new Colour*[_height];
		for (int j = 0; j < _height; j++)
			buffer[i][j] = new Colour[_planes];
	}
	setBoundingBox();
}

Volume::~Volume(void)
{
	clearBuffer();
}

void Volume::clearBuffer(){
	if (buffer != NULL){
		for (int i = 0; i < _width; i++){
			if (buffer[i] != NULL){
				for (int j = 0; j < _height; j++){
					if (buffer[i][j] != NULL){
						//for (int k = 0; k < _planes; k++){
							//if (buffer[i][j][k] != NULL){
							//	delete buffer[i][j][k];
							//}
						//}
						delete[] buffer[i][j];
					}
				}
				delete[] buffer[i];
			}
		}
	}
}

void Volume::allocateBuffer(){
	clearBuffer();
	buffer = new Colour**[_width];
	for (int i = 0; i < _width; i++){
		buffer[i] = new Colour*[_height];
		for (int j = 0; j < _height; j++)
			buffer[i][j] = new Colour[_planes];
	}
}
void Volume::setBoundingBox(){
	/*for (int i = 0; i < _width; i++){
		for (int j = 0; j < _height; j++){
			for (int k = 0; k < _planes; k++){
				if (!buffer[i][j][k].isZero()){
					if (i < boundingBox[0])
						boundingBox[0] = i;
					if (i > boundingBox[1])
						boundingBox[1] = i;
					if (j < boundingBox[0])
						boundingBox[2] = j;
					if (j > boundingBox[1])
						boundingBox[3] = j;
					if (k < boundingBox[0])
						boundingBox[4] = k;
					if (k > boundingBox[1])
						boundingBox[5] = k;
				}
			}
		}
	}*/
	boundingBox[0]= 0;
	boundingBox[2]= 0;
	boundingBox[4]= 0;
	boundingBox[1]= _width;
	boundingBox[3]= _height;
	boundingBox[5]= _planes;
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
	allocateBuffer();

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
	delete[] colourCh;
	delete[] colourF;
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
	

	_width = width;
	_height = height;
	_planes = planes;
	_channels = channels;
	allocateBuffer();

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

	setBoundingBox();
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

bool Volume::hit(const Ray& r, float tmin, float tmax, float time, HitRecord& hit) const{

	return false;
}
