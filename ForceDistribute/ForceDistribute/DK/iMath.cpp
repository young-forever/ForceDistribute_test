#include <vector>
#include "iMath.h"

using namespace std;

int sign(float x)
{
	if (x < 0) return -1;
	else return 1;
}

void  m_cross3d(double *a, double *b, double *c)
{
	// c = axb;  a,b,c all 3 demension;
	double x1 = a[0];
	double y1 = a[1];
	double z1 = a[2];
	double x2 = b[0];
	double y2 = b[1];
	double z2 = b[2];

	c[0] = y1*z2 - y2*z1;
	c[1] = -(x1*z2 - x2*z1);
	c[2] = x1*y2 - x2*y1;
}

vector<vector<double>> matrix_multiply(vector<vector<double>> arrA, vector<vector<double>> arrB)
{

	int rowA = arrA.size();

	int colA = arrA[0].size();

	int rowB = arrB.size();

	int colB = arrB[0].size();

	vector<vector<double>>  res;
	if (colA != rowB)
	{
		return res;
	}
	else
	{
		res.resize(rowA);
		for (int i = 0; i < rowA; ++i)
		{
			res[i].resize(colB);
		}

		for (int i = 0; i < rowA; ++i)
		{
			for (int j = 0; j < colB; ++j)
			{
				for (int k = 0; k < colA; ++k)
				{
					res[i][j] += arrA[i][k] * arrB[k][j];
				}
			}
		}
	}
	return res;
}
