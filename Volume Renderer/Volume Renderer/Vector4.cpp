// Vector4.cpp: implementation of the Vector4 and Matrix4x4 classes.
//
//////////////////////////////////////////////////////////////////////

//#ifdef WIN32
//#include "stdafx.h"
//#else
//#endif
#include "Vector4.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Vector4
//////////////////////////////////////////////////////////////////////

Vector4::Vector4()
{
	v[0]=v[1]=v[2]=v[3]=0;
}

Vector4::Vector4(float x,float y,float z, float w)
{
	v[0]=x;
	v[1]=y;
	v[2]=z;
	v[3]=w;
}

Vector4::Vector4(float *f)
{
	v[0]=f[0];
	v[1]=f[1];
	v[2]=f[2];
	v[3]=f[3];
}

Vector4::Vector4(Vector4D vec)
{
	double *t=vec.VectorTab();
	for (int i=0;i<4;i++)
		v[i]=(float) t[i];
}

Vector4::~Vector4()
{

}


int Vector4::operator==(const Vector4 &vect)
{
	return((v[0]==vect.v[0])&&(v[1]==vect.v[1])&&(v[2]==vect.v[2])&&(v[3]==vect.v[3]));
}

void Vector4::operator=(const Vector4 &vect)
{
	int i;
	for (i=0;i<4;i++) v[i]=vect.v[i];
}

Vector4& Vector4::operator*=(const Matrix4x4 &mat)
{
	Vector4 vect;
	//int i;
	float *m= (float *) &(mat.m[0][0]);
	vect.v[0]=m[0]*v[0]+m[1]*v[1]+m[2]*v[2]+m[3]*v[3];
	vect.v[1]=m[4]*v[0]+m[5]*v[1]+m[6]*v[2]+m[7]*v[3];
	vect.v[2]=m[8]*v[0]+m[9]*v[1]+m[10]*v[2]+m[11]*v[3];
	vect.v[3]=m[12]*v[0]+m[13]*v[1]+m[14]*v[2]+m[15]*v[3];
	(*this)=vect;
	return *this;
}

Vector4 operator*=(float d, const Vector4 &vect)
{
	Vector4 v=vect;
	v*=d;
	return v;
}


Vector4& Vector4::operator+=(const Vector4 &vect)
{
	v[0]+=vect.v[0];
	v[1]+=vect.v[1];
	v[2]+=vect.v[2];
	v[3]+=vect.v[3];
	return (*this);
}

Vector4& Vector4::operator-=(const Vector4 &vect)
{
	v[0]-=vect.v[0];
	v[1]-=vect.v[1];
	v[2]-=vect.v[2];
	v[3]-=vect.v[3];
	return (*this);
}
Vector4 Vector4::operator+(const Vector4 &vect)
{
	Vector4 v=*this;
	v+=vect;
	return (v);
}

Vector4 Vector4::operator*(const float &f)
{
	Vector4 v=*this;
	v[0]*=f;
	v[1]*=f;
	v[2]*=f;
	v[3]*=f;
	return (v);
}

Vector4 Vector4::operator/(const float &f)
{
	Vector4 v=*this;
	float invf=1.f/f;
	v[0]*=invf;
	v[1]*=invf;
	v[2]*=invf;
	v[3]*=invf;
	return (v);
}

Vector4 Vector4::operator-(const Vector4 &vect)
{
	Vector4 v=*this;
	v-=vect;
	return (v);
}

Vector4 Vector4::operator-()
{
	Vector4 v;
	v[0]=-this->v[0];
	v[1]=-this->v[1];
	v[2]=-this->v[2];
	v[3]=-this->v[3];
	return (v);
}

Vector4& Vector4::operator*=(const float &x)
{
	v[0]*=x;
	v[1]*=x;
	v[2]*=x;
	v[3]*=x;
	return (*this);
}


float Vector4::operator*(const Vector4 &vect)
{
	return (v[0]*vect.v[0]+v[1]*vect.v[1]+v[2]*vect.v[2]+v[3]*vect.v[3]);
}

void Vector4::Expand(const Vector4& vect)
{
	v[0]*=vect.v[0];
	v[1]*=vect.v[1];
	v[2]*=vect.v[2];
	v[3]*=vect.v[3];
}

int Vector4::IsWithinSegment(const Vector4 &p1, const Vector4 &p2)
{
Vector4 v1,v2;
v1=p1;
v1-=*this;
v2=p2;
v2-=*this;
//     v1.v[2]=v2.v[2]=v1.v[3]=v2.v[3]=0;
return ((v1*v2)<=0);
}



Vector4 Vector4::CrossProduct3(const Vector4 &v2)
{
	Vector4 result;
	result[0]=(v[1]*v2.v[2]-v[2]*v2.v[1]);
	result[1]=(v[2]*v2.v[0]-v[0]*v2.v[2]);
	result[2]=(v[0]*v2.v[1]-v[1]*v2.v[0]);
	return (result);
}

Vector4 Vector4::operator^(const Vector4 &v2)
{
	Vector4 result;
	result[0]=(v[1]*v2.v[2]-v[2]*v2.v[1]);
	result[1]=(v[2]*v2.v[0]-v[0]*v2.v[2]);
	result[2]=(v[0]*v2.v[1]-v[1]*v2.v[0]);
	return (result);
}

Vector4 Vector4::PlaneEquationFromPoints(const Vector4 &p1,const Vector4 &p2,const Vector4 &p3)
{
	Vector4 v1,v2;
	v1=p2;
	v1-=p1;
	v2=p3;
	v2-=p1;
	v1=v1.CrossProduct3(v2);
	v1.Normalize();
	v1[3]=0;
	v1[3]=-(v1*p1);
	return v1;
}

float Vector4::RayPlaneIntersection(Vector4 p1, Vector4 p2, Vector4 plane)
{
	float t1=p1*plane;
	float t2=p2*plane;
	if (t2-t1!=0.)
		return (-t1)/(t2-t1);
	else return 0.;
}


float *Vector4::VectorTab()
{
	return v;
}

void Vector4::Print(char *s)
{
	printf("%sVecteur 4D: %f %f %f %f \n",s,v[0],v[1],v[2],v[3]);
}

void Vector4::Normalize()
{
	float f;
	f=v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3];
	if (f!=0) f=1/sqrt(f);
	v[0]*=f;v[1]*=f;v[2]*=f;v[3]*=f;
}


float Vector4::Length()
{
	float f;
	f=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
	return f;
}

float Vector4::LengthSquare()
{
	float f;
	f=(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
	return f;
}


Vector4 Vector4::LineProjectionEquation(Vector4 &p0,Vector4 &p1)
{
	Vector4 v,p;
	float l2;

	v=p1-p0;
	v[3]=0; 
	l2=v.LengthSquare();
	if (l2==0.)
		return v;
	l2=1.f/l2;
	v*=l2;
	p=v;
	l2=p0*v;
	p[3]-=l2;
	return p;
}

int Vector4::SkewTest(const Vector4 &line_a,const Vector4 &line_b,
							 const Vector4 &edge_a,const Vector4 &edge_b)
{
	Vector4 v1,v2;
	v1=edge_a;
	v1-=line_a;
	v2=edge_b;
	v2-=line_a;
	v1=v1.CrossProduct3(v2);
	v2=line_b;
	v2-=line_a;	
	return (v1*v2>=0);
}


//////////////////////////////////////////////////////////////////////
//	Matrix4x4
//////////////////////////////////////////////////////////////////////


Matrix4x4::Matrix4x4()
{
	LoadIdentity();
}

Matrix4x4::Matrix4x4(Matrix4x4D mat)
{
	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++)
			m[i][j]=(float) mat.m[i][j];
}

Matrix4x4::~Matrix4x4()
{
	
}

void Matrix4x4::Print(char *s)
{
	printf("%s",s);
	printf("%f %f %f %f\n",  m[0][0], m[0][1], m[0][2], m[0][3]);
	printf("%f %f %f %f\n",  m[1][0], m[1][1], m[1][2], m[1][3]);
	printf("%f %f %f %f\n",  m[2][0], m[2][1], m[2][2], m[2][3]);
	printf("%f %f %f %f\n",  m[3][0], m[3][1], m[3][2], m[3][3]);
}

void Matrix4x4::Inverse()
{

	float m00 = m[0][0];  
	float m01 = m[0][1];  
	float m02 = m[0][2];  
	float m03 = m[0][3];
	float m10 = m[1][0];  
	float m11 = m[1][1];  
	float m12 = m[1][2];  
	float m13 = m[1][3];
	float m20 = m[2][0];  
	float m21 = m[2][1];  
	float m22 = m[2][2];  
	float m23 = m[2][3];
	float m30 = m[3][0];  
	float m31 = m[3][1];  
	float m32 = m[3][2];  
	float m33 = m[3][3];

	float d00 = m11*m22*m33 + m12*m23*m31 + m13*m21*m32 - m31*m22*m13 - m32*m23*m11 - m33*m21*m12;
	float d01 = m10*m22*m33 + m12*m23*m30 + m13*m20*m32 - m30*m22*m13 - m32*m23*m10 - m33*m20*m12;
	float d02 = m10*m21*m33 + m11*m23*m30 + m13*m20*m31 - m30*m21*m13 - m31*m23*m10 - m33*m20*m11;
	float d03 = m10*m21*m32 + m11*m22*m30 + m12*m20*m31 - m30*m21*m12 - m31*m22*m10 - m32*m20*m11;

	float d10 = m01*m22*m33 + m02*m23*m31 + m03*m21*m32 - m31*m22*m03 - m32*m23*m01 - m33*m21*m02;
	float d11 = m00*m22*m33 + m02*m23*m30 + m03*m20*m32 - m30*m22*m03 - m32*m23*m00 - m33*m20*m02;
	float d12 = m00*m21*m33 + m01*m23*m30 + m03*m20*m31 - m30*m21*m03 - m31*m23*m00 - m33*m20*m01;
	float d13 = m00*m21*m32 + m01*m22*m30 + m02*m20*m31 - m30*m21*m02 - m31*m22*m00 - m32*m20*m01;

	float d20 = m01*m12*m33 + m02*m13*m31 + m03*m11*m32 - m31*m12*m03 - m32*m13*m01 - m33*m11*m02;
	float d21 = m00*m12*m33 + m02*m13*m30 + m03*m10*m32 - m30*m12*m03 - m32*m13*m00 - m33*m10*m02;
	float d22 = m00*m11*m33 + m01*m13*m30 + m03*m10*m31 - m30*m11*m03 - m31*m13*m00 - m33*m10*m01;
	float d23 = m00*m11*m32 + m01*m12*m30 + m02*m10*m31 - m30*m11*m02 - m31*m12*m00 - m32*m10*m01;

	float d30 = m01*m12*m23 + m02*m13*m21 + m03*m11*m22 - m21*m12*m03 - m22*m13*m01 - m23*m11*m02;
	float d31 = m00*m12*m23 + m02*m13*m20 + m03*m10*m22 - m20*m12*m03 - m22*m13*m00 - m23*m10*m02;
	float d32 = m00*m11*m23 + m01*m13*m20 + m03*m10*m21 - m20*m11*m03 - m21*m13*m00 - m23*m10*m01;
	float d33 = m00*m11*m22 + m01*m12*m20 + m02*m10*m21 - m20*m11*m02 - m21*m12*m00 - m22*m10*m01;

	float D = m00*d00 - m01*d01 + m02*d02 - m03*d03;
	if (D<0) D=-D;

	if ( D>EPSILON )
	{
		D=1.f/D;
		m[0][0] =  d00*D; 
		m[0][1] = -d10*D;  
		m[0][2] =  d20*D; 
		m[0][3] = -d30*D;
		m[1][0] = -d01*D; 
		m[1][1] =  d11*D;  
		m[1][2] = -d21*D; 
		m[1][3] =  d31*D;
		m[2][0] =  d02*D; 
		m[2][1] = -d12*D;  
		m[2][2] =  d22*D; 
		m[2][3] = -d32*D;
		m[3][0] = -d03*D; 
		m[3][1] =  d13*D;  
		m[3][2] = -d23*D; 
		m[3][3] =  d33*D;
	}
	else
	{
//		std::string text = "\nImpossible to invert a matrice, skipped...\n";
//		throw MathError(IMPOSSIBLE_MATRIX_INVERSION,text,__FILE__,__LINE__);
	}
}

float *Matrix4x4::MatrixTab()
{
	return (float *) m;
}


void Matrix4x4::LoadIdentity()
{
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			m[i][j]=0;
	for (i=0;i<4;i++) m[i][i]=1;
}


void Matrix4x4::operator=(const Matrix4x4 &mat)
{
	int i;
	for (i=0;i<16;i++) m[0][i]=mat.m[0][i];
}

int Matrix4x4::operator!=(const Matrix4x4 &mat)
{
	int i,j,cmp;
	double tmp;
	cmp=1;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
		{
			tmp=m[i][j]-mat.m[i][j];
			if (ABS(tmp)<EPSILON) cmp=0;
		}

	return cmp;
}

Matrix4x4& Matrix4x4::operator*=(const Matrix4x4 &mat)
{
	int i,j,k;
	Matrix4x4 m1;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
		{
			m1.m[i][j]=0;
			for (k=0;k<4;k++)
				m1.m[i][j]+=m[i][k]*mat.m[k][j];
		}
	*this=m1;
	return *this;
}

Matrix4x4& Matrix4x4::operator*=(const float &s)
{
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			m[i][j]*=s;
	return *this;
}

Matrix4x4& Matrix4x4::operator+=(const Matrix4x4 &mat)
{
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			m[i][j]+=mat.m[i][j];
	return *this;

}

Matrix4x4& Matrix4x4::operator-=(const Matrix4x4 &mat)
{
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			m[i][j]-=mat.m[i][j];
	return *this;

}

Matrix4x4 Matrix4x4::operator+(const Matrix4x4 &mat2)
{	
	Matrix4x4 m=*this;
	m+=mat2;
	return m;
}

Matrix4x4 Matrix4x4::operator-(const Matrix4x4 &mat2)
{	
	Matrix4x4 m=*this;
	m-=mat2;
	return m;
}

Matrix4x4 Matrix4x4::operator*(const Matrix4x4 &mat2)
{	
	Matrix4x4 m=*this;
	m*=mat2;
	return m;
}

Matrix4x4 Matrix4x4::operator*(const float &d)
{	
	Matrix4x4 m=*this;
	m*=d;
	return m;
}

Vector4 Matrix4x4::operator*(const Vector4 &v)
{
	Vector4 res;
	res=v;
	res*=*this;
	return res;
}


Matrix4x4 Matrix4x4::operator/(const Matrix4x4 &mat2)
{	
	Matrix4x4 m=*this;
	Matrix4x4 m2=mat2;
	m2.Inverse();
	m*=m2;
	return m;
}

Matrix4x4 operator/(float d, Matrix4x4 &mat2)
{
	Matrix4x4 m=mat2;
	m.Inverse();
	m*=d;
	return m;
}

void Matrix4x4::ZRotate(float angle)
{
	Matrix4x4 tmp;
	tmp.m[0][0]=cos(angle);
	tmp.m[0][1]=-sin(angle);
	tmp.m[1][0]=sin(angle);
	tmp.m[1][1]=cos(angle);
	(*this)*=tmp;
}

void Matrix4x4::YRotate(float angle)
{
	Matrix4x4 tmp;
	tmp.LoadIdentity();
	tmp.m[0][0]=cos(angle);
	tmp.m[0][2]=sin(angle);
	tmp.m[2][0]=-sin(angle);
	tmp.m[2][2]=cos(angle);
	(*this)*=tmp;
}

void Matrix4x4::XRotate(float angle)
{
	Matrix4x4 tmp;
	tmp.LoadIdentity();
	tmp.m[1][1]=cos(angle);
	tmp.m[1][2]=-sin(angle);
	tmp.m[2][1]=sin(angle);
	tmp.m[2][2]=cos(angle);
	(*this)*=tmp;
}
 
void Matrix4x4::Rotate(float angle, Vector4 v)
{
	Rotate(angle,v[0],v[1],v[2]);
}


void Matrix4x4::Rotate(float angle, float x, float y, float z)
{
	Vector4 u,v;
	Matrix4x4 mat,uut,tmp;
	u.Set(x,y,z,0);
	u.Normalize();
	mat[0]=mat[5]=mat[10]=mat[15]=0;
	uut=mat;
	mat[1]=-u[2];mat[2]=u[1];mat[4]=u[2];mat[6]=-u[0];mat[8]=-u[1];mat[9]=u[0];
	uut[0]=u[0]*u[0];uut[1]=u[0]*u[1];uut[2]=u[0]*u[2];
	uut[4]=u[1]*u[0];uut[5]=u[1]*u[1];uut[6]=u[1]*u[2];
	uut[8]=u[2]*u[0];uut[9]=u[2]*u[1];uut[10]=u[2]*u[2];
	tmp-=uut;
	tmp*=cos(angle);
	mat*=sin(angle);
	uut+=tmp;
	uut+=mat;
	uut[15]=1.;
	(*this)*=uut;
}

void Matrix4x4::LookAt(float eyeX ,float  eyeY ,float  eyeZ ,float  centerX ,
				float  centerY ,float  centerZ , 
				float upX ,float  upY ,float  upZ)//might be bugged...
{
	Vector4 f(centerX-eyeX,centerY-eyeY,centerZ-eyeZ,0.);
	f.Normalize();
	Vector4 up(upX,upY,upZ,0.);
	Vector4 upbis=up;
	upbis.Normalize();
	Vector4 s=f.CrossProduct3(upbis);
	Vector4 u=s.CrossProduct3(f);
	Matrix4x4 m;
	m[0]=s[0];m[1]=s[1];m[2]=s[2];m[3]=0;
	m[4]=u[0];m[5]=u[1];m[6]=u[2];m[7]=0;
	m[8]=-f[0];m[9]=-f[1];m[10]=-f[2];m[11]=0;
	(*this)*=m;
	(*this)[12]-=eyeX;
	(*this)[13]-=eyeY;
	(*this)[14]-=eyeZ;
}

void Matrix4x4::Zoom(float coef)
{
	Matrix4x4 m;
	m.LoadIdentity();
	m.m[0][0]=m.m[1][1]=m.m[2][2]=m.m[3][3]=coef;
	*this*=m;
}

void Matrix4x4::Translate(const Vector4 &vect)
{
	Matrix4x4 m;
	m.m[0][3]=vect.v[0];
	m.m[1][3]=vect.v[1];
	m.m[2][3]=vect.v[2];
	*this*=m;
/*
	m[0][4]+=vect.v[0];
	m[1][4]+=vect.v[1];
	m[2][4]+=vect.v[2];
	m[3][4]+=vect.v[3];*/
}

void Matrix4x4::Translate(float dx, float dy, float dz)
{
	Matrix4x4 m;
	m.m[0][3]=dx;
	m.m[1][3]=dy;
	m.m[2][3]=dz;
	*this*=m;
}

void Matrix4x4::Expand(const Vector4 &vect)
{
	Matrix4x4 m;
	m.m[0][0]=vect.v[0];
	m.m[1][1]=vect.v[1];
	m.m[2][2]=vect.v[2];
	*this*=m;
}

void Matrix4x4::Expand(const float dx, const float dy, const float dz)	
{
	Matrix4x4 m;
	m.m[0][0]=dx;
	m.m[1][1]=dy;
	m.m[2][2]=dz;
	(*this)*=m;
}

void Matrix4x4::Tranpose()
{
	int i,j;
	float tmp;
	for (i=0;i<4;i++)
		for (j=0;j<i;j++)
		{
			tmp=m[i][j];
			m[i][j]=m[j][i];
			m[j][i]=tmp;
		}
}

int Matrix4x4::ViewSector(ProjectionPlane p)
{
	switch (p)
	{
	case XY_PLANE :
		return (m[2][0]<0)+(m[2][1]<0)*2+(m[2][2]<0)*4;
	case YZ_PLANE :
		return (m[0][0]<0)+(m[0][1]<0)*2+(m[0][2]<0)*4;
	case ZX_PLANE :
		return (m[1][0]<0)+(m[1][1]<0)*2+(m[1][2]<0)*4;
	default: return -1;
	}
}

void Matrix4x4::Save(FILE *f)
{
	fwrite(m,1,16*sizeof(float),f);
}

void Matrix4x4::Load(FILE *f)
{
	fread(m,1,16*sizeof(float),f);
}

void Matrix4x4::Frustum(float left,float right,float bottom,float top,float near2,float far2)
{
	Matrix4x4 m;
	float a,b,c,d;
	a=(right+left)/(right-left);
	b=(top+bottom)/(top-bottom);
	c=(far2+near2)/(far2-near2);
	d=2.f*far2*near2/(far2-near2);
	m[0]=2*near2/(right-left);
	m[2]=a;
	m[5]=2*near2/(top-bottom);
	m[6]=b;
	m[10]=c;
	m[11]=d;
	m[14]=-1;
	m[15]=0;
	*this*=m;
}

void Matrix4x4::Viewport(int x, int y, int width, int height)
{
	Matrix4x4 m;
	m[0]=width/2.f;
	m[3]=width/2.f+x;
	m[5]=height/2.f;
	m[7]=height/2.f+y;
	*this*=m;
}

void Matrix4x4::MatrixDecomposition(float &rx, float &ry,float &rz,float &tx, float &ty,float &tz)
{
	Matrix4x4 Test;
	double Cx,Cy,Cz,Sx,Sy,Sz;
	Sy=-m[2][0];
	ry=(float) asin(Sy);
	Cy=(float) cos(ry);
	Cx=m[2][2]/Cy;
	Sx=m[2][1]/Cy;
	Cz=m[0][0]/Cy;
	Sz=m[1][0]/Cy;
	if (Cz>0)
		rz=(float) atan(Sz/Cz);
	else
		if (Sz<0)
			rz=(float) (-PI+atan(Sz/Cz));
		else
			rz=(float) (PI+atan(Sz/Cz));
			if (Cx>0)
				rx=(float) atan(Sx/Cx);
			else
				if (Sx<0)
					rx=(float) (-PI+atan(Sx/Cx));
				else
					rx=(float) (PI+atan(Sx/Cx));
	tx=m[0][3];
	ty=m[1][3];
	tz=m[2][3];
	Test.ZRotate(rz);
	Test.YRotate(rx);
	Test.XRotate(ry);

	if ((*this)!=Test)
	{
		if (Sy>0)
			ry=(float) (PI-asin(Sy));
		else
			ry=(float) (-PI-asin(Sy));
		Cy=(float) cos(ry);
		Cx=m[2][2]/Cy;
		Sx=m[2][1]/Cy;
		Cz=m[0][0]/Cy;
		Sz=m[1][0]/Cy;
		if (Cz>0)
			rz=(float) atan(Sz/Cz);
		else
			if (Sz<0)
				rz=(float) (-PI+atan(Sz/Cz));
			else
				rz=(float) (PI+atan(Sz/Cz));
				if (Cx>0)
					rx=(float) atan(Sx/Cx);
				else
					if (Sx<0)
						rx=(float) (-PI-atan(Sx/Cx));
					else
						rx=(float) (PI+atan(Sx/Cx));
	}
}



#undef PI