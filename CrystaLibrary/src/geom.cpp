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
#include "geom.h"


double roundToNthDecimal(double i, int n) {
	return std::round(i * n * 10)/(n * 10);
}

Eigen::Matrix3d getTransformationMatrixFromFractionalToCartesian(std::vector<double> params) {
	
	assert(params.size() == 6);

	Eigen::Matrix3d vector_matrix;
	double a_norm = params[0],
		b_norm = params[1],
		c_norm = params[2],
		alpha = params[3],
		beta = params[4],
		gamma = params[5],
		rad_alpha = alpha * boost::math::constants::pi<double>() / 180,
		rad_beta = beta * boost::math::constants::pi<double>() / 180,
		rad_gamma = gamma * boost::math::constants::pi<double>() / 180;

	double factor, zx_sqrt, Volume;
	//vector a
	vector_matrix.col(0) << a_norm, 0, 0;
	//vector b
	//vector_matrix.col(1) << roundToNthDecimal(b_norm * cos(rad_gamma), ROUND_D), roundToNthDecimal(b_norm * sin(rad_gamma), ROUND_D), 0;
	vector_matrix.col(1) << b_norm * cos(rad_gamma), b_norm * sin(rad_gamma), 0;

	factor = (cos(rad_alpha) - cos(rad_gamma)*cos(rad_beta)) / sin(rad_gamma);
	//zx_sqrt = 1.0 - pow(cos(rad_beta), 2) - pow(factor, 2);
	// zx_sqrt must be higher or equal to 0
	//vector_matrix.col(2) << roundToNthDecimal(c_norm * cos(rad_beta), ROUND_D), roundToNthDecimal(c_norm * factor, ROUND_D), roundToNthDecimal(c_norm * sqrt(abs(zx_sqrt)), ROUND_D);
	//vector_matrix.col(2) << c_norm * cos(rad_beta), c_norm* factor, c_norm* sqrt(abs(zx_sqrt));
	Volume = a_norm * b_norm * c_norm * std::sqrt(1-std::pow(cos(rad_alpha),2)-std::pow(cos(rad_beta), 2)-std::pow(cos(rad_gamma), 2)+ 2 * cos(rad_alpha) * cos(rad_beta) * cos(rad_gamma));
	vector_matrix.col(2) << c_norm * cos(rad_beta), c_norm * factor, Volume / (a_norm*b_norm*sin(rad_gamma));

	return vector_matrix;
}

double getEuclideanDistance(Eigen::VectorXd p1, Eigen::VectorXd p2) {
	double distance = 0;
	int p1_size = p1.size(), p2_size = p2.size();
	assert(p1_size == p2_size);

	for (int i = 0; i < p1_size; ++i) {
		distance += pow(p1(i) - p2(i), 2);
	}
	distance = sqrt(distance);
	return distance;
}

Eigen::Matrix3d getRotationMatrix(Axis a, double angle_degree) {

	double t = angle_degree * M_PI / 180.0, x, y, z;
	Eigen::Matrix3d M;
	// Axis vector
	Eigen::Vector3d N;
	switch (a) {
		case Axis::X :
			N << 1, 0, 0;
			break;
		case Axis::Y:
			N << 0, 1, 0;
			break;
		case Axis::Z:
			N << 0, 0, 1;
			break;
		default:
			break;
	}
	x = N[0]; y = N[1]; z = N[2];
	M.col(0) <<	cos(t) + (1 - cos(t))*pow(x, 2),
		(1 - cos(t))*x*y + z * sin(t),
		(1 - cos(t))*x*z - y*sin(t);

	M.col(1) << (1 - cos(t))*x*y - z * sin(t),
		cos(t) + (1 - cos(t))*pow(y, 2),
		(1 - cos(t))*y*z + x*sin(t);

	M.col(2) << (1 - cos(t))*x*z + y * sin(t),
		(1 - cos(t))*y*z - x * sin(t),
		cos(t) + (1 - cos(t))*pow(z, 2);

	return M;
}

Eigen::Matrix3d getRotationMatrix(Eigen::Vector3d a, double angle_rad) {

	double x, y, z, t = angle_rad;
	Eigen::Matrix3d M;
	// Axis vector
	
	x = a[0]; y = a[1]; z = a[2];
	M.col(0) << cos(t) + (1 - cos(t)) * pow(x, 2),
		(1 - cos(t)) * x * y + z * sin(t),
		(1 - cos(t)) * x * z - y * sin(t);

	M.col(1) << (1 - cos(t)) * x * y - z * sin(t),
		cos(t) + (1 - cos(t)) * pow(y, 2),
		(1 - cos(t)) * y * z + x * sin(t);

	M.col(2) << (1 - cos(t)) * x * z + y * sin(t),
		(1 - cos(t)) * y * z - x * sin(t),
		cos(t) + (1 - cos(t)) * pow(z, 2);

	return M;
}

Eigen::Matrix3d getRotationMatricesComposition(double theta, double phi, double psi) {
	double rad_theta = theta * M_PI / 180.0,
		rad_phi = phi * M_PI / 180.0,
		rad_psi = psi * M_PI / 180.0;
	Eigen::Matrix3d M(3,3);

	M.col(0) << cos(rad_psi)*cos(rad_phi) - sin(rad_psi)*cos(rad_theta)*sin(rad_phi),
		sin(rad_psi)* cos(rad_phi) + cos(rad_psi) * cos(rad_theta) * sin(rad_phi),
		sin(rad_theta) * sin(rad_phi);

	M.col(1) << -cos(rad_psi) * sin(rad_phi) - sin(rad_psi) * cos(rad_theta) * cos(rad_phi),
		-sin(rad_psi)* sin(rad_phi) + cos(rad_psi) * cos(rad_theta) * cos(rad_phi),
		sin(rad_theta)* cos(rad_phi);

	M.col(2) << sin(rad_psi) * sin(rad_theta),
		-cos(rad_psi) * sin(rad_theta),
		cos(rad_theta);

	return M;
}

void processNewFractionalCoordinates(Eigen::VectorXd fractional, Eigen::VectorXd& new_fractional) {

	double k = 0;
	for (int j = 0; j < new_fractional.size(); ++j) {
		// integer part of c(j): the old component
		k = 0;

		// Pick the integer part of the original point
		if (fractional(j) > 0.f) {
			k = floor(fractional(j));
		}
		else {
			if (fractional(j) < 0.f) {
				k = ceil(fractional(j));
			}
		}

		// Processing new fractional coordinates
		// If it is outside: subtract: integer part - original integer part
		// Ex: 1.30 original -> 3.30 new, fix new: 3.30 - (3 - 1) = 1.30
		// Ex (unchanged): 1.30 original -> 1.30 new, fix new: 1.30 - (1 - 1) = 1.30
		/*
		if (new_fractional(j) > 1.f) {
			new_fractional(j) -= (floor(new_fractional(j) - k));
		}
		*/
		// Put inside everything outside of unit cell
		if (new_fractional(j) > 1) {
			new_fractional(j) -= floor(new_fractional(j));
		}
		// If it is negative and outside (could be lower than -1), just add the absolute(floor) to consider the copy inside
		else {
			if (new_fractional(j) < 0) {
				new_fractional(j) += abs(floor(new_fractional(j)));
				//new_fractional(j) = 1.f + new_fractional(j);
			}
		}
	}

}


