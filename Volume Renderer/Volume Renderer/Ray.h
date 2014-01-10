#pragma once

#include "Vector3.h"

class Ray
{
public:
	Ray(void){};
	~Ray(void);
	Ray(const Vector3& a, const Vector3& b){ data[0] = a; data[1] = b;}
	Vector3 origin() const{return data[0];}
	Vector3 direction() const {return data[1];}
	Vector3 pointAtParameter(float t) const {
		return data[0] + t*data[1];
	}

	Vector3 data[2];
};

