#include "forcedistribute.h"
#include<Eigen/Dense>
//#include "mycommand.h"
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
//#include "mainwindow.h"
#include <iostream>
#include "iMath.h"
#include <vector>
#include "main.h"
#include <string>

using namespace std;
using namespace Eigen;
using namespace alglib;

void ExpEndForceCal(VectorXf d_p, VectorXf dd_p, MatrixXf K, MatrixXf B, double *EndForce)//calculate end force.
{

	VectorXf endforce(6);

	endforce = K*d_p + B*dd_p;

	for (short i = 0; i<6; i++)
	{
		EndForce[i] = endforce(i);
		//        cout << EndForce[i]<< endl;
	}
}

string minqplcCal()
{
	//** calculate dir of cables force.
	vector<vector<double>> Rx = { { 1, 0, 0, 0 }, \
	{0, cos(RealPose[3]), -sin(RealPose[3]), 0}, \
	{0, sin(RealPose[3]), cos(RealPose[3]), 0}, \
	{0, 0, 0, 1} };

	vector<vector<double>> Ry = { { cos(RealPose[4]), 0, sin(RealPose[4]), 0 }, \
	{0, 1, 0, 0}, \
	{-sin(RealPose[4]), 0, cos(RealPose[4]), 0}, \
	{0, 0, 0, 1} };

	vector<vector<double>> Rz = { { cos(RealPose[5]), -sin(RealPose[5]), 0, 0 }, \
	{sin(RealPose[5]), cos(RealPose[5]), 0, 0}, \
	{0, 0, 1, 0}, \
	{0, 0, 0, 1} };

	vector<vector<double>> Ptrf = { { 1, 0, 0, RealPose[0] }, { 0, 1, 0, RealPose[1] }, { 0, 0, 1, RealPose[2] }, { 0, 0, 0, 1 } };

	vector<vector<double>> T = matrix_multiply(Ptrf, matrix_multiply(Rz, matrix_multiply(Ry, Rx)));

	vector<vector<double> > V_ep(4, vector<double>(8));

	for (int i = 0; i<4; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			V_ep[i][j] = ep[i][j];
		}
	}
	vector<vector<double>> G_ep = matrix_multiply(T, V_ep);

	double bp1[3] = { bp[0][0], bp[1][0], bp[2][0] };
	double bp2[3] = { bp[0][1], bp[1][1], bp[2][1] };
	double bp3[3] = { bp[0][2], bp[1][2], bp[2][2] };
	double bp4[3] = { bp[0][3], bp[1][3], bp[2][3] };
	double bp5[3] = { bp[0][4], bp[1][4], bp[2][4] };
	double bp6[3] = { bp[0][5], bp[1][5], bp[2][5] };
	double bp7[3] = { bp[0][6], bp[1][6], bp[2][6] };
	double bp8[3] = { bp[0][7], bp[1][7], bp[2][7] };

	double G_ep1[3] = { G_ep[0][0], G_ep[1][0], G_ep[2][0] };
	double G_ep2[3] = { G_ep[0][1], G_ep[1][1], G_ep[2][1] };
	double G_ep3[3] = { G_ep[0][2], G_ep[1][2], G_ep[2][2] };
	double G_ep4[3] = { G_ep[0][3], G_ep[1][3], G_ep[2][3] };
	double G_ep5[3] = { G_ep[0][4], G_ep[1][4], G_ep[2][4] };
	double G_ep6[3] = { G_ep[0][5], G_ep[1][5], G_ep[2][5] };
	double G_ep7[3] = { G_ep[0][6], G_ep[1][6], G_ep[2][6] };
	double G_ep8[3] = { G_ep[0][7], G_ep[1][7], G_ep[2][7] };

	double dir_1[3] = { bp1[0] - G_ep1[0], bp1[1] - G_ep1[1], bp1[2] - G_ep1[2] };
	Map<VectorXd> l_dir_1(dir_1, 3);
	double dir_2[3] = { bp2[0] - G_ep2[0], bp2[1] - G_ep2[1], bp2[2] - G_ep2[2] };
	Map<VectorXd> l_dir_2(dir_2, 3);
	double dir_3[3] = { bp3[0] - G_ep3[0], bp3[1] - G_ep3[1], bp3[2] - G_ep3[2] };
	Map<VectorXd> l_dir_3(dir_3, 3);
	double dir_4[3] = { bp4[0] - G_ep4[0], bp4[1] - G_ep4[1], bp4[2] - G_ep4[2] };
	Map<VectorXd> l_dir_4(dir_4, 3);
	double dir_5[3] = { bp5[0] - G_ep5[0], bp5[1] - G_ep5[1], bp5[2] - G_ep5[2] };
	Map<VectorXd> l_dir_5(dir_5, 3);
	double dir_6[3] = { bp6[0] - G_ep6[0], bp6[1] - G_ep6[1], bp6[2] - G_ep6[2] };
	Map<VectorXd> l_dir_6(dir_6, 3);
	double dir_7[3] = { bp7[0] - G_ep7[0], bp7[1] - G_ep7[1], bp7[2] - G_ep7[2] };
	Map<VectorXd> l_dir_7(dir_7, 3);
	double dir_8[3] = { bp8[0] - G_ep8[0], bp8[1] - G_ep8[1], bp8[2] - G_ep8[2] };
	Map<VectorXd> l_dir_8(dir_8, 3);

	for (int i = 0; i<3; i++)
	{
		dir_1[i] = dir_1[i] / l_dir_1.norm();
		dir_2[i] = dir_2[i] / l_dir_2.norm();
		dir_3[i] = dir_3[i] / l_dir_3.norm();
		dir_4[i] = dir_4[i] / l_dir_4.norm();
		dir_5[i] = dir_5[i] / l_dir_5.norm();
		dir_6[i] = dir_6[i] / l_dir_6.norm();
		dir_7[i] = dir_7[i] / l_dir_7.norm();
		dir_8[i] = dir_8[i] / l_dir_8.norm();
	}

	double r_1[3] = { G_ep1[0] - RealPose[0], G_ep1[1] - RealPose[1], G_ep1[2] - RealPose[2] };
	double r_2[3] = { G_ep2[0] - RealPose[0], G_ep2[1] - RealPose[1], G_ep2[2] - RealPose[2] };
	double r_3[3] = { G_ep3[0] - RealPose[0], G_ep3[1] - RealPose[1], G_ep3[2] - RealPose[2] };
	double r_4[3] = { G_ep4[0] - RealPose[0], G_ep4[1] - RealPose[1], G_ep4[2] - RealPose[2] };
	double r_5[3] = { G_ep5[0] - RealPose[0], G_ep5[1] - RealPose[1], G_ep5[2] - RealPose[2] };
	double r_6[3] = { G_ep6[0] - RealPose[0], G_ep6[1] - RealPose[1], G_ep6[2] - RealPose[2] };
	double r_7[3] = { G_ep7[0] - RealPose[0], G_ep7[1] - RealPose[1], G_ep7[2] - RealPose[2] };
	double r_8[3] = { G_ep8[0] - RealPose[0], G_ep8[1] - RealPose[1], G_ep8[2] - RealPose[2] };

	string x0_s = "[";
	x0_s.append("[");
	x0_s.append(to_string(dir_1[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_2[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_3[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_4[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_5[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_6[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_7[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_8[0]));
	x0_s.append(",");
	x0_s.append(to_string(ExpectEndForce[0]));
	x0_s.append("],");

	x0_s.append("[");
	x0_s.append(to_string(dir_1[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_2[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_3[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_4[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_5[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_6[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_7[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_8[1]));
	x0_s.append(",");
	x0_s.append(to_string(ExpectEndForce[1]));
	x0_s.append("],");

	x0_s.append("[");
	x0_s.append(to_string(dir_1[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_2[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_3[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_4[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_5[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_6[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_7[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_8[2]));
	x0_s.append(",");
	x0_s.append(to_string(ExpectEndForce[2]));
	x0_s.append("],");

	x0_s.append("[");
	x0_s.append(to_string(dir_1[2]*r_1[1] - dir_1[1]*r_1[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_2[2]*r_2[1] - dir_2[1]*r_2[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_3[2]*r_3[1] - dir_3[1]*r_3[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_4[2]*r_4[1] - dir_4[1]*r_4[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_5[2]*r_5[1] - dir_5[1]*r_5[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_6[2]*r_6[1] - dir_6[1]*r_6[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_7[2]*r_7[1] - dir_7[1]*r_7[2]));
	x0_s.append(",");
	x0_s.append(to_string(dir_8[2]*r_8[1] - dir_8[1]*r_8[2]));
	x0_s.append(",");
	x0_s.append(to_string(ExpectEndForce[3]));
	x0_s.append("],");

	x0_s.append("[");
	x0_s.append(to_string(dir_1[0]*r_1[2] - dir_1[2]*r_1[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_2[0]*r_2[2] - dir_2[2]*r_2[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_3[0]*r_3[2] - dir_3[2]*r_3[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_4[0]*r_4[2] - dir_4[2]*r_4[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_5[0]*r_5[2] - dir_5[2]*r_5[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_6[0]*r_6[2] - dir_6[2]*r_6[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_7[0]*r_7[2] - dir_7[2]*r_7[0]));
	x0_s.append(",");
	x0_s.append(to_string(dir_8[0]*r_8[2] - dir_8[2]*r_8[0]));
	x0_s.append(",");
	x0_s.append(to_string(ExpectEndForce[4]));
	x0_s.append("],");

	x0_s.append("[");
	x0_s.append(to_string(dir_1[1]*r_1[0] - dir_1[0]*r_1[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_2[1]*r_2[0] - dir_2[0]*r_2[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_3[1]*r_3[0] - dir_3[0]*r_3[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_4[1]*r_4[0] - dir_4[0]*r_4[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_5[1]*r_5[0] - dir_5[0]*r_5[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_6[1]*r_6[0] - dir_6[0]*r_6[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_7[1]*r_7[0] - dir_7[0]*r_7[1]));
	x0_s.append(",");
	x0_s.append(to_string(dir_8[1]*r_8[0] - dir_8[0]*r_8[1]));
	x0_s.append(",");
	x0_s.append(to_string(ExpectEndForce[5]));
	x0_s.append("]]");

	return x0_s;
}

void EndForceDistrToJointCal(double *expecttesion)//calculate joint tesion.
{
	//
	// This example demonstrates minimization of F(x0,x1) = x0^2 + x1^2 -6*x0 - 4*x1
	// subject to linear constraint x0+x1<=2
	//
	// Exact solution is [x0,x1] = [1.5,0.5]
	//
	// IMPORTANT: this solver minimizes  following  function:
	//     f(x) = 0.5*x'*A*x + b'*x.
	// Note that quadratic term has 0.5 before it. So if you want to minimize
	// quadratic function, you should rewrite it in such way that quadratic term
	// is multiplied by 0.5 too.
	// For example, our function is f(x)=x0^2+x1^2+..., but we rewrite it as 
	//     f(x) = 0.5*(2*x0^2+2*x1^2) + ....
	// and pass diag(2,2) as quadratic term - NOT diag(1,1)!
	//
	real_2d_array a = "[[2,0,0,0,0,0,0,0],[0,2,0,0,0,0,0,0],[0,0,2,0,0,0,0,0],[0,0,0,2,0,0,0,0],[0,0,0,0,2,0,0,0],[0,0,0,0,0,2,0,0],[0,0,0,0,0,0,2,0],[0,0,0,0,0,0,0,2]]";
	real_1d_array b = "[0,0,0,0,0,0,0,0]";
	real_1d_array s = "[1,1,1,1,1,1,1,1]";
	
	string aa = minqplcCal();
	real_2d_array c = aa.c_str();
	integer_1d_array ct = "[0,0,0,0,0,0]";//6个等式约束 类型设置
	real_1d_array x;
	minqpstate state;
	minqpreport rep;

	string x1_s = "[";
	x1_s.append(to_string(ExpectTesion[0]));
	x1_s.append(",");
	x1_s.append(to_string(ExpectTesion[1]));
	x1_s.append(",");
	x1_s.append(to_string(ExpectTesion[2]));
	x1_s.append(",");
	x1_s.append(to_string(ExpectTesion[3]));
	x1_s.append(",");
	x1_s.append(to_string(ExpectTesion[4]));
	x1_s.append(",");
	x1_s.append(to_string(ExpectTesion[5]));
	x1_s.append(",");
	x1_s.append(to_string(ExpectTesion[6]));
	x1_s.append(",");
	x1_s.append(to_string(ExpectTesion[7]));
	x1_s.append("]");

	const char *x0_;
	x0_ = x1_s.c_str();
	real_1d_array x0 = x0_;
	//real_1d_array x0 = "[0,0,0,0,0,0,0,0,0]";

	string x2_s = "[";
	x2_s.append(to_string(MinCalbleForce));
	x2_s.append(",");
	x2_s.append(to_string(MinCalbleForce));
	x2_s.append(",");
	x2_s.append(to_string(MinCalbleForce));
	x2_s.append(",");
	x2_s.append(to_string(MinCalbleForce));
	x2_s.append(",");
	x2_s.append(to_string(MinCalbleForce));
	x2_s.append(",");
	x2_s.append(to_string(MinCalbleForce));
	x2_s.append(",");
	x2_s.append(to_string(MinCalbleForce));
	x2_s.append(",");
	x2_s.append(to_string(MinCalbleForce));
	x2_s.append("]");

	const char *bndl_;
	bndl_ = x2_s.c_str();

	real_1d_array bndl = bndl_;
	//real_1d_array bndl = "[3,3,3,3,3,3,3,3]";
	real_1d_array bndu = "[1000,1000,1000,1000,1000,1000,1000,1000]";

	// create solver, set quadratic/linear terms
	minqpcreate(8, state);
	minqpsetquadraticterm(state, a);
	minqpsetlinearterm(state, b);
	minqpsetlc(state, c, ct);
	minqpsetstartingpoint(state, x0);
	minqpsetbc(state, bndl, bndu);

	// Set scale of the parameters.
	// It is strongly recommended that you set scale of your variables.
	// Knowing their scales is essential for evaluation of stopping criteria
	// and for preconditioning of the algorithm steps.
	// You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
	//
	// NOTE: for convex problems you may try using minqpsetscaleautodiag()
	//       which automatically determines variable scales.
	minqpsetscale(state, s);

	//
	// Solve problem with DENSE-AUL solver.
	//
	// This solver is optimized for problems with up to several thousands of
	// variables and large amount of general linear constraints. Problems with
	// less than 50 general linear constraints can be efficiently solved with
	// BLEIC, problems with box-only constraints can be solved with QuickQP.
	// However, DENSE-AUL will work in any (including unconstrained) case.
	//
	// Default stopping criteria are used.
	//

	minqpsetalgodenseaul(state, 1.0e-5, 1.0e+4, 5);
	minqpoptimize(state);
	minqpresults(state, x, rep);

	//double calTesion[8];
	for (short i = 0; i<8; i++)
	{
		expecttesion[i] = x.operator ()(i);
	}
}

void JointTesionToEndCal(double *realtesion, double *realendforce)
{
	vector<vector<double>> Rx = { { 1, 0, 0, 0 }, \
	{0, cos(RealPose[3]), -sin(RealPose[3]), 0}, \
	{0, sin(RealPose[3]), cos(RealPose[3]), 0}, \
	{0, 0, 0, 1} };

	vector<vector<double>> Ry = { { cos(RealPose[4]), 0, sin(RealPose[4]), 0 }, \
	{0, 1, 0, 0}, \
	{-sin(RealPose[4]), 0, cos(RealPose[4]), 0}, \
	{0, 0, 0, 1} };

	vector<vector<double>> Rz = { { cos(RealPose[5]), -sin(RealPose[5]), 0, 0 }, \
	{sin(RealPose[5]), cos(RealPose[5]), 0, 0}, \
	{0, 0, 1, 0}, \
	{0, 0, 0, 1} };

	vector<vector<double>> Ptrf = { { 1, 0, 0, RealPose[0] }, { 0, 1, 0, RealPose[1] }, { 0, 0, 1, RealPose[2] }, { 0, 0, 0, 1 } };

	vector<vector<double>> T = matrix_multiply(Ptrf, matrix_multiply(Rz, matrix_multiply(Ry, Rx)));

	vector<vector<double> > V_ep(4, vector<double>(8));

	for (int i = 0; i<4; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			V_ep[i][j] = ep[i][j];
		}
	}
	vector<vector<double>> G_ep = matrix_multiply(T, V_ep);

	double bp1[3] = { bp[0][0], bp[1][0], bp[2][0] };
	double bp2[3] = { bp[0][1], bp[1][1], bp[2][1] };
	double bp3[3] = { bp[0][2], bp[1][2], bp[2][2] };
	double bp4[3] = { bp[0][3], bp[1][3], bp[2][3] };
	double bp5[3] = { bp[0][4], bp[1][4], bp[2][4] };
	double bp6[3] = { bp[0][5], bp[1][5], bp[2][5] };
	double bp7[3] = { bp[0][6], bp[1][6], bp[2][6] };
	double bp8[3] = { bp[0][7], bp[1][7], bp[2][7] };

	double G_ep1[3] = { G_ep[0][0], G_ep[1][0], G_ep[2][0] };
	double G_ep2[3] = { G_ep[0][1], G_ep[1][1], G_ep[2][1] };
	double G_ep3[3] = { G_ep[0][2], G_ep[1][2], G_ep[2][2] };
	double G_ep4[3] = { G_ep[0][3], G_ep[1][3], G_ep[2][3] };
	double G_ep5[3] = { G_ep[0][4], G_ep[1][4], G_ep[2][4] };
	double G_ep6[3] = { G_ep[0][5], G_ep[1][5], G_ep[2][5] };
	double G_ep7[3] = { G_ep[0][6], G_ep[1][6], G_ep[2][6] };
	double G_ep8[3] = { G_ep[0][7], G_ep[1][7], G_ep[2][7] };

	double dir_1[3] = { bp1[0] - G_ep1[0], bp1[1] - G_ep1[1], bp1[2] - G_ep1[2] };
	Map<VectorXd> l_dir_1(dir_1, 3);
	double dir_2[3] = { bp2[0] - G_ep2[0], bp2[1] - G_ep2[1], bp2[2] - G_ep2[2] };
	Map<VectorXd> l_dir_2(dir_2, 3);
	double dir_3[3] = { bp3[0] - G_ep3[0], bp3[1] - G_ep3[1], bp3[2] - G_ep3[2] };
	Map<VectorXd> l_dir_3(dir_3, 3);
	double dir_4[3] = { bp4[0] - G_ep4[0], bp4[1] - G_ep4[1], bp4[2] - G_ep4[2] };
	Map<VectorXd> l_dir_4(dir_4, 3);
	double dir_5[3] = { bp5[0] - G_ep5[0], bp5[1] - G_ep5[1], bp5[2] - G_ep5[2] };
	Map<VectorXd> l_dir_5(dir_5, 3);
	double dir_6[3] = { bp6[0] - G_ep6[0], bp6[1] - G_ep6[1], bp6[2] - G_ep6[2] };
	Map<VectorXd> l_dir_6(dir_6, 3);
	double dir_7[3] = { bp7[0] - G_ep7[0], bp7[1] - G_ep7[1], bp7[2] - G_ep7[2] };
	Map<VectorXd> l_dir_7(dir_7, 3);
	double dir_8[3] = { bp8[0] - G_ep8[0], bp8[1] - G_ep8[1], bp8[2] - G_ep8[2] };
	Map<VectorXd> l_dir_8(dir_8, 3);

	for (int i = 0; i<3; i++)
	{
		dir_1[i] = dir_1[i] / l_dir_1.norm();
		dir_2[i] = dir_2[i] / l_dir_2.norm();
		dir_3[i] = dir_3[i] / l_dir_3.norm();
		dir_4[i] = dir_4[i] / l_dir_4.norm();
		dir_5[i] = dir_5[i] / l_dir_5.norm();
		dir_6[i] = dir_6[i] / l_dir_6.norm();
		dir_7[i] = dir_7[i] / l_dir_7.norm();
		dir_8[i] = dir_8[i] / l_dir_8.norm();
	}
	realendforce[0] = realtesion[0] * dir_1[0] + realtesion[1] * dir_2[0] + realtesion[2] * dir_3[0] + realtesion[3] * dir_4[0] + \
		realtesion[4] * dir_5[0] + realtesion[5] * dir_6[0] + realtesion[6] * dir_7[0] + realtesion[7] * dir_8[0];
	realendforce[1] = realtesion[0] * dir_1[1] + realtesion[1] * dir_2[1] + realtesion[2] * dir_3[1] + realtesion[3] * dir_4[1] + \
		realtesion[4] * dir_5[1] + realtesion[5] * dir_6[1] + realtesion[6] * dir_7[1] + realtesion[7] * dir_8[1];
	realendforce[2] = realtesion[0] * dir_1[2] + realtesion[1] * dir_2[2] + realtesion[2] * dir_3[2] + realtesion[3] * dir_4[2] + \
		realtesion[4] * dir_5[2] + realtesion[5] * dir_6[2] + realtesion[6] * dir_7[2] + realtesion[7] * dir_8[2];

	double r_1[3] = { G_ep1[0] - RealPose[0], G_ep1[1] - RealPose[1], G_ep1[2] - RealPose[2] };
	double r_2[3] = { G_ep2[0] - RealPose[0], G_ep2[1] - RealPose[1], G_ep2[2] - RealPose[2] };
	double r_3[3] = { G_ep3[0] - RealPose[0], G_ep3[1] - RealPose[1], G_ep3[2] - RealPose[2] };
	double r_4[3] = { G_ep4[0] - RealPose[0], G_ep4[1] - RealPose[1], G_ep4[2] - RealPose[2] };
	double r_5[3] = { G_ep5[0] - RealPose[0], G_ep5[1] - RealPose[1], G_ep5[2] - RealPose[2] };
	double r_6[3] = { G_ep6[0] - RealPose[0], G_ep6[1] - RealPose[1], G_ep6[2] - RealPose[2] };
	double r_7[3] = { G_ep7[0] - RealPose[0], G_ep7[1] - RealPose[1], G_ep7[2] - RealPose[2] };
	double r_8[3] = { G_ep8[0] - RealPose[0], G_ep8[1] - RealPose[1], G_ep8[2] - RealPose[2] };

	double Tor1[3], Tor2[3], Tor3[3], Tor4[3], Tor5[3], Tor6[3], Tor7[3], Tor8[3];
	double Fc_1[3], Fc_2[3], Fc_3[3], Fc_4[3], Fc_5[3], Fc_6[3], Fc_7[3], Fc_8[3];

	Fc_1[0] = realtesion[0] * dir_1[0];
	Fc_1[1] = realtesion[0] * dir_1[1];
	Fc_1[2] = realtesion[0] * dir_1[2];
	Fc_2[0] = realtesion[1] * dir_2[0];
	Fc_2[1] = realtesion[1] * dir_2[1];
	Fc_2[2] = realtesion[1] * dir_2[2];
	Fc_3[0] = realtesion[2] * dir_3[0];
	Fc_3[1] = realtesion[2] * dir_3[1];
	Fc_3[2] = realtesion[2] * dir_3[2];
	Fc_4[0] = realtesion[3] * dir_4[0];
	Fc_4[1] = realtesion[3] * dir_4[1];
	Fc_4[2] = realtesion[3] * dir_4[2];
	Fc_5[0] = realtesion[4] * dir_5[0];
	Fc_5[1] = realtesion[4] * dir_5[1];
	Fc_5[2] = realtesion[4] * dir_5[2];
	Fc_6[0] = realtesion[5] * dir_6[0];
	Fc_6[1] = realtesion[5] * dir_6[1];
	Fc_6[2] = realtesion[5] * dir_6[2];
	Fc_7[0] = realtesion[6] * dir_7[0];
	Fc_7[1] = realtesion[6] * dir_7[1];
	Fc_7[2] = realtesion[6] * dir_7[2];
	Fc_8[0] = realtesion[7] * dir_8[0];
	Fc_8[1] = realtesion[7] * dir_8[1];
	Fc_8[2] = realtesion[7] * dir_8[2];

	m_cross3d(r_1, Fc_1, Tor1);
	m_cross3d(r_2, Fc_2, Tor2);
	m_cross3d(r_3, Fc_3, Tor3);
	m_cross3d(r_4, Fc_4, Tor4);
	m_cross3d(r_5, Fc_5, Tor5);
	m_cross3d(r_6, Fc_6, Tor6);
	m_cross3d(r_7, Fc_7, Tor7);
	m_cross3d(r_8, Fc_8, Tor8);

	realendforce[3] = Tor1[0] + Tor2[0] + Tor3[0] + Tor4[0] + Tor5[0] + Tor6[0] + Tor7[0] + Tor8[0];
	realendforce[4] = Tor1[1] + Tor2[1] + Tor3[1] + Tor4[1] + Tor5[1] + Tor6[1] + Tor7[1] + Tor8[1];
	realendforce[5] = Tor1[2] + Tor2[2] + Tor3[2] + Tor4[2] + Tor5[2] + Tor6[2] + Tor7[2] + Tor8[2];
}
