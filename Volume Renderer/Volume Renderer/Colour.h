#pragma once

#include <ostream>

#define RGB 3
#define RGBA 4
#define GREYSCALE 1
#define COMPRGB 6
#define COMPGREY 2

class Colour
{
public:

	Colour(void)
	{
		_nbChannels = 3;
		_intensities = new float[3];
		for (int i = 0; i < _nbChannels; i++){
			_intensities[i] = 0;
		}
	}

	Colour(const int nbChannels, const float intensities[])
	{
		_nbChannels = nbChannels;
		_intensities = new float[nbChannels];
		for (int i = 0; i < nbChannels; i++){
			_intensities[i] = intensities[i];
		}
	}

	Colour(const Colour& original) {*this = original;}

	~Colour(void)
	{
		delete _intensities;
	}

	bool isRGB() const{return _nbChannels == RGB;}
	bool isRGBA() const{return _nbChannels == RGBA;};
	bool isGreyscale() const{return _nbChannels == GREYSCALE;};
	bool isCOMPRGB() const{return _nbChannels == COMPRGB;};
	bool isCOMPGREY() const{return _nbChannels == COMPGREY;};
	bool isZero() const;
	
	int nbChannels() const {return _nbChannels;}
	float intensity(int c) const {return _intensities[c];}
	float operator[](int c){return _intensities[c];}
	void setIntensity(int c, float i) {_intensities[c] = i;}

	Colour& operator=(const Colour & right_op);
	Colour& operator+=(const Colour & right_op);
	Colour& operator*=(const Colour & right_op);
	Colour& operator/=(const Colour & right_op);
	Colour& operator*(float right_op);
	Colour& operator/(float right_op);

	Colour operator+() const{return *this;}
	Colour operator-()const{
		Colour newC(_nbChannels, _intensities);
		for (int i = 0; i < _nbChannels; i++) newC._intensities[i] = -_intensities[i];
		return newC;}

	friend std::ostream& operator<<( std::ostream & out, const Colour & theColour);
	friend Colour operator*(const Colour& c, float f);
	friend Colour operator*(float f, const Colour& c);
	friend Colour operator/(const Colour& c, float f);
	friend Colour operator*(const Colour& c1, const Colour& c2);
	friend Colour operator/(const Colour& c1, const Colour& c2);
	friend Colour operator+(const Colour& c1, const Colour& c2);

	void clamp();

	int _nbChannels;
	float *_intensities;


};

inline Colour& Colour::operator=(const Colour & right_op) {
	_nbChannels = right_op._nbChannels;
	if (_intensities != NULL)
		delete _intensities;
	_intensities = new float[_nbChannels];
	for (int i = 0; i < _nbChannels; i++){
		_intensities[i] = right_op._intensities[i];
	}
	return *this;
}

inline Colour& Colour::operator+=(const Colour & right_op){
	if (_nbChannels != right_op._nbChannels)
		throw "incompatible colours";
	for (int i = 0; i < _nbChannels; i++)
		_intensities[i] += right_op._intensities[i];
	return *this;
}
inline Colour& Colour::operator*=(const Colour & right_op){
	if (_nbChannels != right_op._nbChannels)
		throw "incompatible colours";
	for (int i = 0; i < _nbChannels; i++)
		_intensities[i] *= right_op._intensities[i];
	return *this;
}
inline Colour& Colour::operator/=(const Colour & right_op){
	if (_nbChannels != right_op._nbChannels)
		throw "incompatible colours";
	for (int i = 0; i < _nbChannels; i++)
		_intensities[i] /= right_op._intensities[i];
	return *this;
}
inline Colour& Colour::operator*(float right_op){
	for (int i = 0; i < _nbChannels; i++)
		_intensities[i] *= right_op;
	return *this;
}
inline Colour& Colour::operator/(float right_op){
	for (int i = 0; i < _nbChannels; i++)
		_intensities[i] /= right_op;
	return *this;
}

inline std::ostream& operator<<( std::ostream & out, const Colour & theColour){
	for (int i = 0; i < theColour._nbChannels; i++){
		out << theColour._intensities[i] << ' ';
	}
	out << '/n';
	return out;
}

inline Colour operator*(const Colour& c, float f){
	Colour newC(c._nbChannels, c._intensities);
	for (int i = 0; i < newC._nbChannels; i++) 
		newC._intensities[i] *= f ;
	return newC;
}
inline Colour operator*(float f, const Colour& c){
	Colour newC(c._nbChannels, c._intensities);
	for (int i = 0; i <  newC._nbChannels; i++) 
		newC._intensities[i] *= f ;
	return newC;
}
inline Colour operator/(const Colour& c, float f){
	Colour newC(c._nbChannels, c._intensities);
	for (int i = 0; i <  newC._nbChannels; i++) 
		newC._intensities[i] /= f ;
	return newC;
}
inline Colour operator*(const Colour& c1, const Colour& c2){
	if (c1._nbChannels != c2._nbChannels)
		throw "incompatible colours";
	Colour newC(c1._nbChannels, c1._intensities);
	for (int i = 0; i <  newC._nbChannels; i++) 
		newC._intensities[i] *= c2._intensities[i] ;
	return newC;
}
inline Colour operator/(const Colour& c1, const Colour& c2){
	if (c1._nbChannels != c2._nbChannels)
		throw "incompatible colours";
	Colour newC(c1._nbChannels, c1._intensities);
	for (int i = 0; i < newC._nbChannels; i++) 
		newC._intensities[i] /= c2._intensities[i] ;
	return newC;
}
inline Colour operator+(const Colour& c1, const Colour& c2){
	if (c1._nbChannels != c2._nbChannels)
		throw "incompatible colours";
	Colour newC(c1._nbChannels, c1._intensities);
	for (int i = 0; i <  newC._nbChannels; i++) 
		newC._intensities[i] += c2._intensities[i] ;
	return newC;
}

inline void Colour::clamp(){
	for (int i = 0; i < _nbChannels; i++){
		if (_intensities[i] >1.0f)
			_intensities[i] = 1.0f;
		if (_intensities[i] < 0.0f)
			_intensities[i] = 0.0f;
	}
}

inline bool Colour::isZero()const{
	bool isZero = true;
	for (int k = 0; k < _nbChannels; k++)
		if (_intensities[k] != 0.0f)
			isZero = false;
	return isZero;
}
