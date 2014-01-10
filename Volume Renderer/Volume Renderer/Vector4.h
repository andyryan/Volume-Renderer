// Vector4.h: interface for the Vector4 class.
//
//////////////////////////////////////////////////////////////////////

#ifndef Vector4_H
#define Vector4_H

#define DELTA_FLOAT 0.00001
#define EPSILON 1.E-14f
#define EPSILOND 1.E-34f
//#define FHUGE 1.E+34f
#define FHUGE 0x7f800000
#define DHUGE 1.E+36f
#define PI 3.141592

//#include "DirectTraceDef.h"
#include <stdio.h>
#include "Vector4D.h"

enum ProjectionPlane
{
	XY_PLANE,
	YZ_PLANE,
	ZX_PLANE
};


class Matrix4x4;

class Vector4  
{
private:
	friend class Vector4D;
	friend class Matrix4x4;
public:
	float v[4];
	
	Vector4();
	Vector4(float x,float y,float z,float w);
	Vector4(float *v);
	Vector4(Vector4D v);
	virtual ~Vector4();
	void operator=(const Vector4 &v);
	Vector4& operator*=(const Matrix4x4 &mat);
	Vector4 operator*(const float &f);
	float operator*(const Vector4 &vect);
	Vector4& operator*=(const float &x);
	friend Vector4 operator*(double d, Vector4& mat2);
	Vector4& operator+=(const Vector4 &vect);
	Vector4& operator-=(const Vector4 &vect);
	Vector4 operator+(const Vector4 &vect);
	Vector4 operator-(const Vector4 &vect);
	Vector4 operator-();
	int operator==(const Vector4 &vect);
	Vector4 operator/(const float &f);
	float &operator[](const int &i) {return v[i];};
	void Normalize();
	int IsWithinSegment(const Vector4 &v1, const Vector4 &v2);
	void Expand(const Vector4& vect);
	inline void Set(float x, float y, float z, float w) {v[0]=x;v[1]=y;v[2]=z;v[3]=w;};
	inline void Set(float vect[4]) { v[0]=vect[0];v[1]=vect[1];v[2]=vect[2];v[3]=vect[3];};
	Vector4 CrossProduct3(const Vector4 &v2);
	Vector4 operator^(const Vector4 &v2);
	void Print(char *s);
	float *VectorTab();
	float Length();
	float LengthSquare();
	static Vector4 LineProjectionEquation(Vector4 &p0,Vector4 &p1); 
	static float RayPlaneIntersection(Vector4 p1, Vector4 p2, Vector4 plane);  
	static Vector4 PlaneEquationFromPoints(const Vector4 &p1,const Vector4 &p2,const Vector4 &p3);
	static int SkewTest(const Vector4 &line_a,const Vector4 &line_b,const Vector4 &edge_a, const Vector4 &edge_b);
};

typedef Vector4 Point4;

class Matrix4x4
{
private:
	friend class Vector4;
	friend class Matrix4x4D;
public:
	float m[4][4];
	Matrix4x4();
	Matrix4x4(Matrix4x4D mat);
	~Matrix4x4();
	void operator=(const Matrix4x4 &mat);
	float& operator[](const int &i) {return m[i>>2][i&3];};
	Matrix4x4& operator*=(const Matrix4x4 &mat);
	Matrix4x4& operator*=(const float &s);
	Matrix4x4& operator+=(const Matrix4x4 &mat);
	Matrix4x4& operator-=(const Matrix4x4 &mat);
	Matrix4x4 operator+(const Matrix4x4 &mat2);
	Matrix4x4 operator-(const Matrix4x4 &mat2);
	Matrix4x4 operator*(const Matrix4x4 &mat2);
	Vector4 operator*(const Vector4 &v);
	Matrix4x4 operator*(const float &d);
	Matrix4x4 operator/(const Matrix4x4 &mat2);
	friend Matrix4x4 operator/(float d, Matrix4x4 &mat2);
	int operator!=(const Matrix4x4 &mat);

	void Inverse();
	void ZRotate(float angle);
	void YRotate(float angle);
	void XRotate(float angle);
	void Rotate(float angle, float x, float y, float z);
	void Rotate(float angle, Vector4 v);
	void LookAt(float eyeX ,float  eyeY ,float  eyeZ ,float  centerX ,
				float  centerY ,float  centerZ , 
				float upX ,float  upY ,float  upZ);
	void Translate(const Vector4 &vect);
	void Translate(float dx, float dy, float dz);
	void Tranpose();
	void Expand(const Vector4 &vect);	
	void Expand(const float dx, const float dy, const float dz);	
	void Zoom(float coef);
	void LoadIdentity();
	void MatrixDecomposition(float &rx, float &ry,float &rz,float &tx, float &ty,float &tz);
	void Frustum(float left,float right,float bottom,float top,float near,float far);
	void Viewport(int x, int y, int width, int height);
	int ViewSector(ProjectionPlane p);
	void Print(char *s);
	void Save(FILE *f);
	void Load(FILE *f);
	float *MatrixTab();
};


//typedef PAN_FLOAT Vector3D[3]; 

/*------------------------macros concernant les vecteurs--------------------------------------------*/ 



#define MIN(a,b)      (((a)<(b))?(a):(b))

#define MAX(a,b)      (((a)>(b))?(a):(b))

#define LERP(a,l,h)   ((l)+(((h)-(l))*(a)))

#define Vscale(S,a)   (a)[0]*=S; (a)[1]*=S; (a)[2]*=S;

#define Vnegate(a)    (a)[0]=0.0-(a)[0];\
			(a)[1]=0.0-(a)[1];\
			(a)[2]=0.0-(a)[2]

#define Vdot(a,b)     ((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])

#define Vcross(a,b,c)	 (c)[0]=(a)[1]*(b)[2]-(a)[2]*(b)[1];\
			 (c)[1]=(a)[2]*(b)[0]-(a)[0]*(b)[2];\
			 (c)[2]=(a)[0]*(b)[1]-(a)[1]*(b)[0]

#define Vlen(a)       (sqrt(Vdot(a,a)))

#define Vcopy(a,b)	 (b)[0]=(a)[0];(b)[1]=(a)[1];(b)[2]=(a)[2];

#define Vadd(a,b,c)	 (c)[0]=(a)[0]+(b)[0];\
			 (c)[1]=(a)[1]+(b)[1];\
			 (c)[2]=(a)[2]+(b)[2]

#define Vsub(a,b,c)	 (c)[0]=(a)[0]-(b)[0];\
			 (c)[1]=(a)[1]-(b)[1];\
			 (c)[2]=(a)[2]-(b)[2]

#define Vcomb(A,a,B,b,c)	(c)[0]=(A)*(a)[0]+(B)*(b)[0];\
				(c)[1]=(A)*(a)[1]+(B)*(b)[1];\
			 	(c)[2]=(A)*(a)[2]+(B)*(b)[2]

#define VaddS(A,a,b,c)	 (c)[0]=(A)*(a)[0]+(b)[0];\
				 (c)[1]=(A)*(a)[1]+(b)[1];\
				 (c)[2]=(A)*(a)[2]+(b)[2]

#define Vprint(msg,v)		printf("%s %g %g %g\n", msg,\
					(v)[0],(v)[1],(v)[2])

#define Vzero(v)	(v)[0]=0.0;(v)[1]=0.0;v[2]=0.0

#define Vnull(a) (	( (a)[0] < EPSILON) && ( (a)[0] > -EPSILON) && ( (a)[1] < EPSILON) && ( (a)[1] > -EPSILON) && ( (a)[2] < EPSILON) && ( (a)[2] > -EPSILON)   ) 

#define Varequal(a,b) ( ((a)[0]==(b)[0]) && ((a)[1]==(b)[1]) && ((a)[2]==(b)[2]))

#define Vareclose(a,b) ( (ABS((a)[0]-(b)[0])<=FEPSILON) && (ABS((a)[1]-(b)[1])<=FEPSILON) && (ABS((a)[2]-(b)[2])<=FEPSILON))

#define ABS(a)	   (((a) < 0) ? -(a) : (a))

#endif


