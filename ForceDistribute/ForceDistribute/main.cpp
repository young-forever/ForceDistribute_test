#include <iostream>
#include "forcedistribute.h"
#include "main.h"
using namespace std;

double MinCalbleForce = 5.0;//绳索力下极限；
double RealPose[6] = { 0,0,0.3,0.1,0.1,0.1 };
double bp[4][8] = { { -0.406,	-0.393,	 0.401,	0.387,	-0.405,	-0.395,	 0.400,	0.390 }, 
					{0.289,		-0.309,	-0.295,	0.302,	0.289,	-0.306,	-0.294,	0.303} ,
					{0.670,		0.672,	0.670,	0.670,	0.005,	0.005,	0.005,	0.005},
					{1,1,1,1,1,1,1,1} };
double ep[4][8] = { {-0.042,	-0.042,	0.016,	0.016,	-0.016,	-0.016,	0.042,	0.042},
					{0.031,		-0.031,	-0.050,	0.050,	0.050,	-0.050,	-0.031,	0.031},
					{-0.026,	-0.026,	-0.026,	-0.026,	0.026,	0.026,	0.026,	0.026},
					{1,1,1,1,1,1,1,1} };

double ExpectTesion[8] = { 5,5,5,5,19.03,16.45,17.37,16.145};
double ExpectEndForce[6] = { 0, 0, 5, 0.1, 0, 0.3 };

int main()
{
	EndForceDistrToJointCal(ExpectTesion);//力分配函数，输入末端力状态，输出期望绳索力
	
	double testforce[6];
	JointTesionToEndCal(ExpectTesion, testforce);

	cout << "ExpectTesion: ";
	for (int i = 0; i < 8;i++)
	{
		cout << ExpectTesion[i] << "   ";
	}
	cout << endl;

	cout << "Endforce: ";
	for (int i = 0; i < 6;i++)
	{
		cout << testforce[i] << "   ";
	}
	cout << endl;

	while (1)
		;

	return 0;
}



