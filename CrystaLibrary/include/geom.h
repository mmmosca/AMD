/*
 * Permission is granted to copy, distribute and/or modify the documents
 * in this directory and its subdirectories unless otherwise stated under
 * the terms of the GNU Free Documentation License, Version 1.1 or any later version 
 * published by the Free Software Foundation; with no Invariant Sections, 
 * no Front-Cover Texts and no Back-Cover Texts. A copy of the license 
 * is available at the website of the GNU Project.
 * The programs and code snippets in this directory and its subdirectories
 * are free software; you can redistribute them and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later
 * version. This code is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * Author Marco M. Mosca, email: marcomichele.mosca@gmail.com
*/
#ifndef _GEOM_H
#define _GEOM_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <assert.h>
#include <queue>
#include <map>

#include <boost/math/special_functions.hpp>
#include <boost/math/constants/constants.hpp>

#define M_PI           3.14159265358979323846
#define ROUND_D	20
#define PRINTSTEPINFO(step, r1, v1, r2, v2, r3, v3)	std::cout << "|--> Step " << (step) << ":\n" << \
													"\t" << to_string((r1)) << ":\t" << (v1)->transpose() << ",\n" << \
													"\t" << to_string((r2)) << ":\t" << (v2)->transpose() << ",\n" << \
													"\t" << to_string((r3)) << ":\t" << (v3)->transpose() << std::endl;
#define PRINTREDUCEINFO(message, r1, r2, k)	std::cout << "\t|--> " << message << ": " << (r1) << "->" << (r2) << ", " << (k) << std::endl;
#define math_sign(X) ((X) < 0.f) ? -1 : 1
#define MAX(X,Y)	((X)>=(Y)) ? (X) : (Y)
#define MIN(X,Y)	((X)<=(Y)) ? (X) : (Y)

enum class Axis {X,Y,Z};

using namespace std;

double roundToNthDecimal(double i, int n);

/* 
*	This function returns the Cartesian vector components of the 3 unit cell axis 
*	param1 : vector of lengths (a, b, c) and angles (alpha, beta, gamma)
*/
Eigen::Matrix3d getTransformationMatrixFromFractionalToCartesian(std::vector<double> params);

double getEuclideanDistance(Eigen::VectorXd p1, Eigen::VectorXd p2);

Eigen::Matrix3d getRotationMatrix(Axis a, double angle_degree);

Eigen::Matrix3d getRotationMatrix(Eigen::Vector3d a, double angle_rad );

Eigen::Matrix3d getRotationMatricesComposition(double theta, double phi, double psi);

void processNewFractionalCoordinates(Eigen::VectorXd fractional, Eigen::VectorXd& new_fractional);

/********************\
|**** POLYHEDRON ****|
\********************/


#endif // !_GEOM_H
