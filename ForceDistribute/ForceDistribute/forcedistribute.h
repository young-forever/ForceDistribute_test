#ifndef FORCEDISTRIBUTE_H
#define FORCEDISTRIBUTE_H
#include<Eigen/Dense>
using namespace Eigen;

void ExpEndForceCal(VectorXf d_p, VectorXf dd_p, MatrixXf K, MatrixXf B, double *EndForce);//calculate end force.
void EndForceDistrToJointCal(double *expecttesion);//calculate joint tesion.
void JointTesionToEndCal(double *realtesion, double *realendforce);

#endif // FORCEDISTRIBUTE_H
