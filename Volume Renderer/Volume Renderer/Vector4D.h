// Vector4D.h: interface for the Vector4D class.
//
//////////////////////////////////////////////////////////////////////

#ifndef Vector4D_H
#define Vector4D_H

//#include "DirectTraceDef.h"
#include <stdio.h>
#include "Vector4.h"


class Matrix4x4D;

class Vector4D  
{
private:
	friend class Matrix4x4D;
	friend class Vector4;
public:
	double v[4];
	
	Vector4D();
	Vector4D(double x,double y,double z,double w);
	Vector4D(double *v);
	Vector4D(Vector4 v);
	virtual ~Vector4D();
	void operator=(const Vector4D &v);
	Vector4D& operator*=(const Matrix4x4D &mat);
	Vector4D operator*(const double &f);
	double operator*(const Vector4D &vect);
	Vector4D& operator*=(const double &x);
	friend Vector4D operator*(double d, Vector4D& mat2);
	Vector4D& operator+=(const Vector4D &vect);
	Vector4D& operator-=(const Vector4D &vect);
	Vector4D operator+(const Vector4D &vect);
	Vector4D operator-(const Vector4D &vect);
	Vector4D operator-();
	int operator==(const Vector4D &vect);
	Vector4D operator/(const double &f);
	double &operator[](const int &i) {return v[i];};
	void Normalize();
	int IsWithinSegment(const Vector4D &v1, const Vector4D &v2);
	void Expand(const Vector4D& vect);
	void Set(double x, double y, double z, double w);
	void Set(double vect[4]);
	Vector4D CrossProduct3(const Vector4D &v2);
	Vector4D operator^(const Vector4D &v2);
	void Print(char *s);
	double *VectorTab();
	double Length();
	double LengthSquare();
	static Vector4D LineProjectionEquation(Vector4D &p0,Vector4D &p1); 
	static double RayPlaneIntersection(Vector4D p1, Vector4D p2, Vector4D plane);  
	static Vector4D PlaneEquationFromPoints(const Vector4D &p1,const Vector4D &p2,const Vector4D &p3);
	static int SkewTest(const Vector4D &line_a,const Vector4D &line_b,const Vector4D &edge_a, const Vector4D &edge_b);
};

typedef Vector4D Point4D;

class Matrix4x4D
{
private:
	friend class Vector4D;
	friend class Matrix4x4;
public:
	double m[4][4];
	Matrix4x4D();
	Matrix4x4D(Matrix4x4 mat);
	~Matrix4x4D();

	void operator=(const Matrix4x4D &mat);
	double& operator[](const int &i) {return m[i>>2][i&3];};
	Matrix4x4D& operator*=(const Matrix4x4D &mat);
	Matrix4x4D& operator*=(const double &s);
	Matrix4x4D& operator+=(const Matrix4x4D &mat);
	Matrix4x4D& operator-=(const Matrix4x4D &mat);
	Matrix4x4D operator+(const Matrix4x4D &mat2);
	Matrix4x4D operator-(const Matrix4x4D &mat2);
	Matrix4x4D operator*(const Matrix4x4D &mat2);
	Vector4D operator*(const Vector4D &v);
	Matrix4x4D operator*(const double &d);
	Matrix4x4D operator/(const Matrix4x4D &mat2);
	friend Matrix4x4D operator/(double d, Matrix4x4D &mat2);
	int operator!=(const Matrix4x4D &mat);

	void Inverse();
	void ZRotate(double angle);
	void YRotate(double angle);
	void XRotate(double angle);
	void Rotate(double angle, double x, double y, double z);
	void Rotate(double angle, Vector4D v);
	void LookAt(double eyeX ,double  eyeY ,double  eyeZ ,double  centerX ,
				double  centerY ,double  centerZ , 
				double upX ,double  upY ,double  upZ);
	void Translate(const Vector4D &vect);
	void Translate(double dx, double dy, double dz);
	void Tranpose();
	void Expand(const Vector4D &vect);	
	void Expand(const double dx, const double dy, const double dz);	
	void Zoom(double coef);
	void LoadIdentity();
	void MatrixDecomposition(double &rx, double &ry,double &rz,double &tx, double &ty,double &tz);
	void Frustum(double left,double right,double bottom,double top,double near,double far);
	void Viewport(int x, int y, int width, int height);
	//int ViewSector(ProjectionPlane p);
	void Print(char *s);
	void Save(FILE *f);
	void Load(FILE *f);
	double *MatrixTab();
};



#endif

