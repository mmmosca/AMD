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
#include <periodicpointcloud.h>

/*****************************************\
|********  PERIODIC POINT CLOUD  *********|
\*****************************************/

std::vector<Eigen::VectorXd>& PeriodicPointCloud::getPointCloud() {
	return this->G_combinations;
}

void PeriodicPointCloud::updateCoefficientsForEveryDirection(int n) {
	//assert(n > 0 && "Number of positive coefficients must be > 0");
	this->clearCoefficients();
	for (int i = -n; i <= n; ++i) {
		this->C.push_back(i);
		++this->coef_count;
	}
}

void PeriodicPointCloud::updateCoefficientsForPositiveDirection(int n) {
	//assert(n > 0 && "Number of positive coefficients must be > 0");
	this->clearCoefficients();
	for (int i = 0; i <= n; ++i) {
		this->C.push_back(i);
		++this->coef_count;
	}
}

void PeriodicPointCloud::updateCoefficientsFromRange(int i, int j) {
	//assert(n > 0 && "Number of positive coefficients must be > 0");
	this->clearCoefficients();
	for (int n = i; n <= j; ++n) {
		this->C.push_back(n);
		++this->coef_count;
	}
}

bool checkElementsInVector(Eigen::VectorXd main_vector, std::vector<int> to_check) {
	for (int c : to_check) {
		for (int m = 0; m < main_vector.size(); ++m) {
			if (c == main_vector[m]) return true;
		}
	}
	return false;
}

void PeriodicPointCloud::updateCoefficientCombinations(std::vector<int> to_filter, Eigen::VectorXd v, int n) {
	if (n == this->generator_count) {
		return;
	}
	else {
		if (n == -1) {
			//Inizialize the starting vector and call the recursion with n = 0
			Eigen::VectorXd first_combination(this->generator_count);
			for (int i = 0; i < this->generator_count; ++i) {
				first_combination(i) = this->C[0];
			}
			this->C_combinations.push_back(first_combination);
			++this->c_combination_count;
			updateCoefficientCombinations(to_filter, first_combination, n + 1);

		}
		else {
			updateCoefficientCombinations(to_filter, v, n + 1);
			Eigen::VectorXd next_combination(v);
			for (int c = 1; c < this->coef_count; ++c) {
				next_combination[n] = this->C[c];
				if (to_filter.empty()) {					
					this->C_combinations.push_back(next_combination);
					++this->c_combination_count;
				}
				else {
					if (checkElementsInVector(next_combination, to_filter)) {
						this->C_combinations.push_back(next_combination);
						++this->c_combination_count;
					}
				}
				updateCoefficientCombinations(to_filter, next_combination, n + 1);
			}
		}
	}
}

void PeriodicPointCloud::updateGeneratorCombinations() {
	int i, j;
	Eigen::VectorXd c;

	// If linear combinations should be added start from new ones
	if (this->g_combination_count > 0) {
		j = this->g_combination_count;
	}
	else {
		// If linear combinations should be created from beginning (from first coefficient permutation)
		j = 0;
	}
	for (; j < this->c_combination_count; ++j) {
		c = this->C_combinations[j];
		Eigen::VectorXd new_combination = Eigen::VectorXd::Zero(this->dimension);
		i = 0;
		for (auto &g : this->G) {
			new_combination += (c[i] * g);
			++i;
		}
		this->G_combinations.push_back(new_combination);
		++this->g_combination_count;
	}
}

void PeriodicPointCloud::updateGeneratorCombinationsForPoints(std::vector<Eigen::VectorXd> points) {
	int i, j;
	Eigen::VectorXd c;

	// If linear combinations should be added start from new ones
	if (this->g_combination_count > 0) {
		j = this->g_combination_count;
	}
	else {
		// If linear combinations should be created from beginning (from first coefficient permutation)
		j = 0;
	}

	for (; j < this->c_combination_count; ++j) {
		c = this->C_combinations[j];
		Eigen::VectorXd new_combination = Eigen::VectorXd::Zero(this->dimension);
		i = 0;
		for (auto &g : this->G) {
			new_combination += (c[i] * g);
			++i;
		}
		for (auto &p : points) {
			Eigen::VectorXd new_motif(new_combination + p);
			this->G_combinations.push_back(new_motif);
			++this->g_combination_count;
		}
	}
}

void PeriodicPointCloud::extendBasisVectorsCoefficientsAndLinearCombinations(int n, std::vector<Eigen::VectorXd> points) {
	// Assumption: this->C is ordered
	int left_border = this->C[0], right_border = this->C[this->coef_count - 1];
	std::vector<int> new_coefficients;

	for (int i = left_border - n; i < left_border; ++i) {
		new_coefficients.push_back(i);
	}
	for (int i = right_border + 1; i <= right_border + n; ++i) {
		new_coefficients.push_back(i);
	}

	left_border -= n;
	right_border += n;

	// Update the internal coefficient list with new ones and keep the order
	this->updateCoefficientsFromRange(left_border, right_border);
	// Update coefficient permutations
	this->updateCoefficientCombinations(new_coefficients);
	if (points.empty()) {
		this->updateGeneratorCombinations();
	}
	else {
		this->updateGeneratorCombinationsForPoints(points);
	}
}

void PeriodicPointCloud::addGenerator(Eigen::VectorXd g) {
	assert(g.size() == this->dimension);
	this->G.insert(g);
	++this->generator_count;
}

void PeriodicPointCloud::clear() {
	this->C.clear();
	this->G.clear();
	this->C_combinations.clear();
	this->G_combinations.clear();
	this->generator_count = 0;
	this->coef_count = 0;
	this->c_combination_count = 0;
	this->g_combination_count = 0;
}

void PeriodicPointCloud::clearCoefficients() {
	this->C.clear();
	this->coef_count = 0;
}

void PeriodicPointCloud::clearCombinations() {
	//this->C.clear();
	//this->G.clear();
	this->C_combinations.clear();
	this->G_combinations.clear();
	//this->generator_count = 0;
	//this->coef_count = 0;
	this->c_combination_count = 0;
	this->g_combination_count = 0;

}

void PeriodicPointCloud::print_info() {
	std::cout << "Generators: " << this->generator_count << std::endl;
	for (auto gen : this->G) {
		std::cout << '\t' << gen.transpose() << std::endl;
	}
	std::cout << "Coefficients: " << this->coef_count << std::endl;
	std::cout << '\t';
	for (auto c : this->C) {
		std::cout << c << " ";
	}
	std::cout << std::endl;
	std::cout << "Coefficients combinations: " << this->c_combination_count << std::endl;
	for (auto elem : this->C_combinations) {
		std::cout << '\t' << elem.transpose() << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Generators combinations: " << this->g_combination_count << std::endl;
	for (auto elem : this->G_combinations) {
		std::cout << '\t' << elem.transpose() << std::endl;
	}

	return;
}
