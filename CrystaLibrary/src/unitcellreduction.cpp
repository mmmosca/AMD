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
#include <unitcellreduction.h>


bool reduceVector_Buerger(Eigen::Vector3d* toReduce, Eigen::Vector3d* refVector, Eigen::Vector3d** shortest, int *matrix_entry) {
	bool reduced = false;
	Eigen::Vector3d ref = *refVector;
	int size = ref.size(), sign;
	double cos_redref, norm_toReduce, norm_ref, fractional_proj;

	norm_toReduce = toReduce->norm();
	norm_ref = ref.norm();
	cos_redref = toReduce->dot(ref) / (norm_toReduce * norm_ref);
	sign = (cos_redref < 0) ? 1 : -1;
	fractional_proj = roundToNthDecimal(abs((norm_toReduce / norm_ref) * cos_redref), 10);

	while (fractional_proj > 0.5) {
		for (int i = 0; i < size; ++i) {
			(*toReduce)(i) = (*toReduce)(i) + sign * ref(i);
		}
		norm_toReduce = toReduce->norm();
		cos_redref = toReduce->dot(ref) / (norm_toReduce * norm_ref);
		fractional_proj = roundToNthDecimal(abs((norm_toReduce / norm_ref) * cos_redref), 10);
		++(*matrix_entry);
		reduced = true;
	}
	(*matrix_entry) *= sign;
	*shortest = (toReduce->norm() <= ref.norm()) ? toReduce : refVector;
	return reduced;
}

bool reduceVector_Niggli_old(Eigen::Vector3d* toReduce, Eigen::Vector3d* refVector, Eigen::Vector3d* waiting, bool& face_diagonal_check, Eigen::Vector3d** shortest, int *matrix_entry) {
	bool reduced = false;
	Eigen::Vector3d ref = *refVector;
	int size = ref.size(), sign;
	double cos_redref, norm_toReduce, norm_ref, fractional_proj, ref_waiting_angle_character, toreduce_waiting_character;

	face_diagonal_check = false;
	norm_toReduce = toReduce->norm();
	norm_ref = ref.norm();
	cos_redref = toReduce->dot(ref) / (norm_toReduce * norm_ref);
	sign = (cos_redref < 0) ? 1 : -1;
	// Round the fractional projection to avoid floating point problem in the while guard
	fractional_proj = roundToNthDecimal((norm_toReduce / norm_ref) * cos_redref, 20);
	ref_waiting_angle_character = 2 * ref.dot(*waiting);
	toreduce_waiting_character = 2 * toReduce->dot(*waiting);
	// Add a condition when the swapped couple have been already tested, if their fractional_proj was 0.5 and this is the viceversa case -> do not reduce
	// In tuples add case 4: it comes the couple of vector have been already reduced and the viceversa should be tested
	while ((abs(fractional_proj) > 0.5) || ((fractional_proj == 0.5) && (2 * toreduce_waiting_character < ref_waiting_angle_character)) || ((fractional_proj == -0.5) && (ref_waiting_angle_character < 0))) {
		if (abs(fractional_proj) == 0.5) {
			// Fix to avoid over-reduction with face-diagonal case. It could continue to inf since abs(fractional_proj) == 0.5 for two consequent iterations.
			if (face_diagonal_check) break;
			else face_diagonal_check = true;
			//std::cout << "Face diagonal case" << std::endl;
		}
		for (int i = 0; i < size; ++i) {
			(*toReduce)(i) = (*toReduce)(i) + sign * ref(i);
		}
		
		norm_toReduce = toReduce->norm();
		cos_redref = toReduce->dot(ref) / (norm_toReduce * norm_ref);
		fractional_proj = roundToNthDecimal((norm_toReduce / norm_ref) * cos_redref, 20);

		++(*matrix_entry);
		reduced = true;
	}
	(*matrix_entry) *= sign;
	*shortest = (toReduce->norm() <= ref.norm()) ? toReduce : refVector;
	return reduced;
}

bool reduceVector_Niggli(Eigen::Vector3d* toReduce, Eigen::Vector3d* refVector, Eigen::Vector3d* waiting, bool& face_diagonal_check, int *matrix_entry) {
	bool reduced = false;
	Eigen::Vector3d ref = *refVector;
	int size = ref.size(), sign;
	double cos_redref, norm_toReduce, norm_ref, fractional_proj, ref_waiting_angle_character, toreduce_waiting_character;

	face_diagonal_check = false;
	norm_toReduce = toReduce->norm();
	norm_ref = ref.norm();
	cos_redref = toReduce->dot(ref) / (norm_toReduce * norm_ref);
	sign = (cos_redref < 0) ? 1 : -1;
	// Round the fractional projection to avoid floating point problem in the while guard
	fractional_proj = roundToNthDecimal((norm_toReduce / norm_ref) * cos_redref, 50);
	//fractional_proj = (norm_toReduce / norm_ref) * cos_redref;
	ref_waiting_angle_character = 2 * ref.dot(*waiting);
	toreduce_waiting_character = 2 * toReduce->dot(*waiting);
	// Add a condition when the swapped couple have been already tested, if their fractional_proj was 0.5 and this is the viceversa case -> do not reduce
	// In tuples add case 4: it comes the couple of vector have been already reduced and the viceversa should be tested
	while ((abs(fractional_proj) > 0.5) || ((fractional_proj == 0.5) && (2 * toreduce_waiting_character < ref_waiting_angle_character)) || ((fractional_proj == -0.5) && (ref_waiting_angle_character < 0))) {
		//std::cout << "INSIDE LOOP" << std::endl;
		if (abs(fractional_proj) == 0.5) {
			// Fix to avoid over-reduction with face-diagonal case. It could continue to inf since abs(fractional_proj) == 0.5 for two consequent iterations.
			if (face_diagonal_check) break;
			else face_diagonal_check = true;
			//std::cout << "Face diagonal case" << std::endl;
		}
		for (int i = 0; i < size; ++i) {
			(*toReduce)(i) = (*toReduce)(i) + sign * ref(i);
		}

		norm_toReduce = toReduce->norm();
		cos_redref = toReduce->dot(ref) / (norm_toReduce * norm_ref);
		fractional_proj = roundToNthDecimal((norm_toReduce / norm_ref) * cos_redref, 50);
		//fractional_proj = (norm_toReduce / norm_ref) * cos_redref;

		++(*matrix_entry);
		reduced = true;
	}
	(*matrix_entry) *= sign;
	return reduced;
}

struct EigenVCompareNorm {
	bool operator() (const Eigen::Vector3d* u, const Eigen::Vector3d* v) { return u->norm() < v->norm(); }
} EigenNormCompare;

bool redirectAxis(std::vector<Eigen::Vector3d*> basis, std::map<Eigen::Vector3d*, int> mapPointerToIndex, Eigen::Matrix3d& transform) {
	Eigen::Vector3d* b1 = basis[0], *b2 = basis[1], *b3 = basis[2];
	bool redirected = false;
	// signs of angle cos (angle character) between vectors
	double cos_b1b2, cos_b1b3, cos_b2b3;
	int sign_cos_b1b2, sign_cos_b1b3, sign_cos_b2b3,
		// entries of the transformation matrix to rediret vectors: 1 -> identity vector component
		b1_entry = 1, b2_entry = 1, b3_entry = 1,
		// multiplication of angles signs to distinguish cases
		orientation,
		// reference to a right angle in the system
		right = -1,
		// integer references to vector basis
		b1_ref, b2_ref, b3_ref;

	b1_ref = mapPointerToIndex[b1];
	b2_ref = mapPointerToIndex[b2];
	b3_ref = mapPointerToIndex[b3];

	// Look at the angle characters (signs)
	//	0->right
	//	1->acute
	//	-1->obtuse
	cos_b1b2 = b1->dot(*b2);
	sign_cos_b1b2 = math_sign(cos_b1b2);
	cos_b1b3 = b1->dot(*b3);
	sign_cos_b1b3 = math_sign(cos_b1b3);
	cos_b2b3 = b2->dot(*b3);
	sign_cos_b2b3 = math_sign(cos_b2b3);

	orientation = sign_cos_b1b2 * sign_cos_b1b3 * sign_cos_b2b3;

	// orientation == 1
	// case 1: 2 obtuse, 1 acute -> all-acute
	// Solution: Look for the vector opposite to the acute angle and change its orientation

	// case 2: all acute -> no change from this if guard	
	if (orientation == 1) {
		// Turn the orientation of b1: obtuse -> acute
		if (cos_b2b3 < 0) {
			redirected = true;
			b1_entry = -1;
#ifdef DEBUG
			PRINTREDUCEINFO("RD", b1_ref, "acute", b1_entry);
#endif
		};
		// Turn the orientation of b2: obtuse -> acute
		if (cos_b1b3 < 0) {
			redirected = true;
			b2_entry = -1;
#ifdef DEBUG
			PRINTREDUCEINFO("RD", b2_ref, "acute", b2_entry);
#endif
		};
		// Turn the orientation of b3: obtuse -> acute
		if (cos_b1b2 < 0) {
			redirected = true;
			b3_entry = -1;
#ifdef DEBUG
			PRINTREDUCEINFO("RD", b3_ref, "acute", b3_entry);
#endif
		};


	}
	// All angles obtuse -> no changes
	else if ((sign_cos_b1b2 == -1) && (sign_cos_b1b3 == -1) && (sign_cos_b2b3 == -1)) {
		// leave all obtuse
	}

	// orientation == -1
	// case 1: 1 obtuse, 2 acute -> all obtuse
	// Solution: Look for the vector opposite to the obtuse angle and change its orientation
	// case 2: all obtuse -> no change

	// orientation == 0
	// case 1: Save the right angle with cos = 0
	else {
		//if ((orientation == 0 || orientation == -1))
		// Turn the orientation of b3: acute -> obtuse
		if (cos_b1b2 > 0) {
			redirected = true;
			b3_entry = -1;
#ifdef DEBUG
			PRINTREDUCEINFO("RD", b3_ref, "obtuse", b3_entry);
#endif
		}
		// b3 opposite angle is 90 degrees: change its orientation
		if (cos_b1b2 == 0) {
			right = b3_ref;
#ifdef DEBUG
			PRINTREDUCEINFO("RD-right", b3_ref, "obtuse", b3_entry);
#endif
		};
		// Turn the orientation of b2: acute -> obtuse
		if (cos_b1b3 > 0) {
			redirected = true;
			b2_entry = -1;
#ifdef DEBUG
			PRINTREDUCEINFO("RD", b2_ref, "obtuse", b2_entry);
#endif
		}
		// b2 opposite angle is 90 degrees: change its orientation
		if (cos_b1b3 == 0) {
			right = b2_ref;
#ifdef DEBUG
			PRINTREDUCEINFO("RD-right", b2_ref, "obtuse", b2_entry);
#endif
		};
		// Turn the orientation of b1: acute -> obtuse
		if (cos_b2b3 > 0) {
			redirected = true;
			b1_entry = -1;
#ifdef DEBUG
			PRINTREDUCEINFO("RD", b1_ref, "obtuse", b1_entry);
#endif
		}
		// b1 opposite angle is 90 degrees: change its orientation
		if (cos_b2b3 == 0) {
			right = b1_ref;
#ifdef DEBUG
			PRINTREDUCEINFO("RD-right", b1_ref, "obtuse", b1_entry);
#endif
		};
	}

	int right_entry = -1;
	// Update the transformation matrix
	transform(b1_ref, b1_ref) = b1_entry;
	transform(b2_ref, b2_ref) = b2_entry;
	transform(b3_ref, b3_ref) = b3_entry;
	// case right angles
	// If (there is at least one right angle) &&
	// (one of them has been changed to obtuse) (possible cases 1 * 1 * -1, 1 * -1 * 1, -1 * 1 * 1) 
	// then redirect the vector opposite to the right angle
	// Solution: with right angles redirect opposite vectors (to get all acute)
	if ((right != -1) && (b1_entry * b2_entry * b3_entry == -1)) { redirected = true; transform(right, right) = right_entry; }

	// Update the vector basis with the previous changes
	for (int i = 0; i < b1->size(); ++i) {
		if (right == b1_ref) {
			(*b1)(i) = right_entry * (*b1)(i);
		}
		else {
			(*b1)(i) = b1_entry * (*b1)(i);
		}
		if (right == b2_ref) {
			(*b2)(i) = right_entry * (*b2)(i);
		}
		else {
			(*b2)(i) = b2_entry * (*b2)(i);
		}
		if (right == b3_ref) {
			(*b3)(i) = right_entry * (*b3)(i);
		}
		else {
			(*b3)(i) = b3_entry * (*b3)(i);
		}
	}
	return redirected;
}

Eigen::Matrix3d reduceUnitCell_old(std::vector<double> &cell_parameters, Eigen::Matrix3d &transform, bool &total_reduced) {

	// matrices for axis of the new system, temporary transformation matrix, final transformation matrix
	Eigen::Matrix3d new_axis, transf_matrix, final_transf_matrix;
	// Old axis system
	Eigen::Matrix3d axis_matrix = getTransformationMatrixFromFractionalToCartesian(cell_parameters);
	std::vector<Eigen::Vector3d*> axis_vector(3);
	// 0: A, 1: B, 2: C
	axis_vector[0] = new Eigen::Vector3d(axis_matrix.col(0));
	axis_vector[1] = new Eigen::Vector3d(axis_matrix.col(1));
	axis_vector[2] = new Eigen::Vector3d(axis_matrix.col(2));
	// Map from axis pointer to its id (as above)
	std::map<Eigen::Vector3d*, int> mapPointerToIndex;
	mapPointerToIndex.insert(std::pair<Eigen::Vector3d*, int>(axis_vector[0], 0));
	mapPointerToIndex.insert(std::pair<Eigen::Vector3d*, int>(axis_vector[1], 1));
	mapPointerToIndex.insert(std::pair<Eigen::Vector3d*, int>(axis_vector[2], 2));
	// queue where vector to be reduce are pushed 
	std::queue<std::tuple<Eigen::Vector3d*, Eigen::Vector3d*, Eigen::Vector3d*, int>> reducing_room;
	// pointer to: vector to be reduced, reference vector, shortest after a reduction, unchanged vector
	Eigen::Vector3d *toReduce, *refVector, *shortest, *waiting;
	// int ids of vectors, matrix entry to transform the vector to be reduced
	int ref, red, sh, wait, matrix_entry;
	bool reduced, face_diagonal;
	//std::cout << "Old axis (columns):\n" << axis_matrix << std::endl;

	final_transf_matrix.setIdentity();
	std::sort(axis_vector.begin(), axis_vector.end(), EigenNormCompare);
	// Start to reduce the longest with respect to the smallest
	reducing_room.push(std::make_tuple(axis_vector[1], axis_vector[0], axis_vector[2], 1));
	//std::cout << "\t|--> Reduce X with respect to (->) Y" << std::endl;
	while (!reducing_room.empty()) {
		//	Get the top element reference
		std::tuple<Eigen::Vector3d*, Eigen::Vector3d*, Eigen::Vector3d*, int> current_tuple = reducing_room.front();
		//	Remove the element
		reducing_room.pop();
		shortest = nullptr;
		matrix_entry = 0;
		// Get pointer to vector of the top element
		toReduce = std::get<0>(current_tuple);
		refVector = std::get<1>(current_tuple);
		waiting = std::get<2>(current_tuple);

		ref = mapPointerToIndex[refVector];
		red = mapPointerToIndex[toReduce];
		wait = mapPointerToIndex[waiting];

		//Redirect the axis to converge to type 1 or type 2 cell (all acute or all obtuse angles)	
		transf_matrix.setIdentity();
		redirectAxis(axis_vector, mapPointerToIndex, transf_matrix);
		final_transf_matrix = final_transf_matrix * transf_matrix;

		// To Do: detect fractional_proj == +-0.5
		// Add parameter true or false
		reduced = reduceVector_Niggli_old(toReduce, refVector, waiting, face_diagonal, &shortest, &matrix_entry);
		//std::cout << face_diagonal << std::endl;
		//reduced = reduceVector(toReduce, refVector, &shortest, &matrix_entry);

		sh = mapPointerToIndex[shortest];

		// Update the transformation matrix
		if (reduced) {
			transf_matrix.setIdentity();
			// Negate the entry to transform the motif
			transf_matrix(red, ref) = -matrix_entry;
			final_transf_matrix = final_transf_matrix * transf_matrix;
		}

		// Save information about cell that is being reduced or not
		total_reduced = total_reduced || reduced;
#ifdef DEBUG
		if (reduced) { PRINTREDUCEINFO("Check Y" , red, ref, -matrix_entry); }
		else { PRINTREDUCEINFO("Check N" , red, ref, ""); }
#endif

		//	SWITCH: Check the pointers of the first and the second and check the integer code:
		//		0: No push -> just check the second condition of a couple
		//		1: check 1st couple of vectors
		//		2: check 2nd couple of vectors
		//		3: check 3rd couple of vectors
		//		4: if fractional_proj == +-0.5 and came from the previous 0 case -> do not push
		switch (std::get<3>(current_tuple)) {
			case 0: {
				// No push -> just check the second condition of a couple
				// face_diagonal avoid pushing a case where reduction loop never stops
				// The opposite couple have been reduced and it was in face diagonal case -> Do not push
				if (reduced) {
					// A vector has been reduced, so remove the next couple -> it is invalid
					if (!reducing_room.empty()) reducing_room.pop();
					reducing_room.push(std::make_tuple(refVector, toReduce, waiting, 1));
				}
				else {
					// Couple fixed -> no push
				}
				break;
			}
			case 1: {
				if (reduced) {
					//PRINTSTEPINFO(1, ref, refVector, red, toReduce, wait, waiting);
					reducing_room.push(std::make_tuple(refVector, toReduce, waiting, 1));
				}
				else {
					// First couple -> Check the viceversa and push the shortest one with a new one
					if (toReduce == shortest) {
						//PRINTSTEPINFO(1, wait, waiting, sh, shortest, ref, refVector);
						reducing_room.push(std::make_tuple(refVector, toReduce, waiting, 0));
						reducing_room.push(std::make_tuple(waiting, shortest, refVector, 2));
					}
					else {
						//PRINTSTEPINFO(1, wait, waiting, sh, shortest, red, toReduce);
						reducing_room.push(std::make_tuple(refVector, toReduce, waiting, 0));
						reducing_room.push(std::make_tuple(waiting, shortest, toReduce, 2));
					}
				}
				break;
			}
			case 2: {
				if (reduced) {
					//PRINTSTEPINFO(2, ref, refVector, red, toReduce, wait, waiting);
					reducing_room.push(std::make_tuple(refVector, toReduce, waiting, 1));
				}
				else {
					// Second couple -> Check the viceversa and push the shortest one with the last one
					//PRINTSTEPINFO(2, red, toReduce, wait, waiting, ref, refVector);
					reducing_room.push(std::make_tuple(refVector, toReduce, waiting, 0));
					reducing_room.push(std::make_tuple(toReduce, waiting, refVector, 3));
				}
				break;
			}
			case 3: {
				// Last step: first 2 couple have been fixed -> check the last combination
				if (reduced) {
					//PRINTSTEPINFO(3, ref, refVector, red, toReduce, wait, waiting);
					// Enqueue to start from the beginning to check if they need to be reduced each other
					reducing_room.push(std::make_tuple(refVector, toReduce, waiting, 1));
				}
				else {
					//std::cout << "\t|--> End" << std::endl;
					// Push just to check the viceversa
					reducing_room.push(std::make_tuple(refVector, toReduce, waiting, 0));
				}

				break;
			}
			default: {
				break;
			}
		}

		// Check the body-diagonal only after the reduction of every vector
		if (reducing_room.empty()) {
			Eigen::Vector3d body_diagonal = Eigen::Vector3d(*(axis_vector[0]) + *(axis_vector[1]) + *(axis_vector[2]));
			std::sort(axis_vector.begin(), axis_vector.end(), EigenNormCompare);
			Eigen::Vector3d b_face_diagonal = *(axis_vector[0]) + *(axis_vector[2]), 
				a_face_diagonal = *(axis_vector[1]) + *(axis_vector[2]);
			// Guard 1: || a + b + c || < || c || or
			// The body-diagonal is smaller then the longest vector -> replace c with the body-diagonal
			// Guard 2: || a + b + c || == || c || and || a + c || > || b + c || 
			// body diagonal is equal to longest vector -> check the face-diagonals with c
			if ( (body_diagonal.norm() < axis_vector[2]->norm()) || 
					( (body_diagonal.norm() == axis_vector[2]->norm()) && (b_face_diagonal.norm() > a_face_diagonal.norm()) ) ) {
				
				// Replace c with the body-diagonal
				std::map<Eigen::Vector3d*, int>::iterator pointer_int_it = mapPointerToIndex.find(axis_vector[2]);
				int axis_id = pointer_int_it->second;
							
				(*axis_vector[2])[0] = body_diagonal[0];
				(*axis_vector[2])[1] = body_diagonal[1];
				(*axis_vector[2])[2] = body_diagonal[2];

				// The 3rd vector is the 3rd + 2nd + 1st (row 2 = 1 1 1) for the axis
				// Points are transformed by row 2 = (-1,-1,1)
				transf_matrix.setIdentity();
				Eigen::Vector3d body_diag_vector = Eigen::Vector3d(-1.f, -1.f, -1.f);
				body_diag_vector[axis_id] = 1.f;
				transf_matrix.row(axis_id) = body_diag_vector;
				final_transf_matrix = final_transf_matrix * transf_matrix;
				std::sort(axis_vector.begin(), axis_vector.end(), EigenNormCompare);
				reducing_room.push(std::make_tuple(axis_vector[1], axis_vector[0], axis_vector[2], 1));
#ifdef DEBUG
				PRINTREDUCEINFO("BD", axis_id, "Body-Diagonal", "");
#endif
			}
		}
	}

	// Update cell parameters and the transformation matrix
	for (auto v : mapPointerToIndex) {
		new_axis.col(v.second) = *(v.first);
	}
	cell_parameters[0] = new_axis.col(0).norm();
	cell_parameters[1] = new_axis.col(1).norm();
	cell_parameters[2] = new_axis.col(2).norm();

	cell_parameters[3] = acos(new_axis.col(1).dot(new_axis.col(2)) / (new_axis.col(1).norm() * new_axis.col(2).norm())) * 180.f / M_PI;
	cell_parameters[4] = acos(new_axis.col(0).dot(new_axis.col(2)) / (new_axis.col(0).norm() * new_axis.col(2).norm())) * 180.f / M_PI;
	cell_parameters[5] = acos(new_axis.col(0).dot(new_axis.col(1)) / (new_axis.col(0).norm() * new_axis.col(1).norm())) * 180.f / M_PI;

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			transform(i, j) = final_transf_matrix(i, j);
		}
	}

	return new_axis;
	//std::cout << "New axis (columns):\n" << new_axis << std::endl;
	//std::cout << "Transformation matrix:\n" << transform << std::endl;
}

Eigen::Matrix3d reduceUnitCell(std::vector<double> &cell_parameters, Eigen::Matrix3d &transform, bool &total_reduced) {

	total_reduced = false;
	// matrices for axis of the new system, temporary transformation matrix, final transformation matrix
	Eigen::Matrix3d new_axis, transf_matrix, final_transf_matrix;
	// Old axis system
	Eigen::Matrix3d axis_matrix = getTransformationMatrixFromFractionalToCartesian(cell_parameters);
	std::vector<Eigen::Vector3d*> axis_vector(3);
	// 0: A, 1: B, 2: C
	axis_vector[0] = new Eigen::Vector3d(axis_matrix.col(0));
	axis_vector[1] = new Eigen::Vector3d(axis_matrix.col(1));
	axis_vector[2] = new Eigen::Vector3d(axis_matrix.col(2));
	// Map from axis pointer to its id (as above)
	std::map<Eigen::Vector3d*, int> mapPointerToIndex;
	mapPointerToIndex.insert(std::pair<Eigen::Vector3d*, int>(axis_vector[0], 0));
	mapPointerToIndex.insert(std::pair<Eigen::Vector3d*, int>(axis_vector[1], 1));
	mapPointerToIndex.insert(std::pair<Eigen::Vector3d*, int>(axis_vector[2], 2));
	// queue where vector to be reduce are pushed 
	std::queue<std::tuple<Eigen::Vector3d*, Eigen::Vector3d*, Eigen::Vector3d*, int>> reducing_room;
	// pointer to: vector to be reduced, reference vector, shortest after a reduction, unchanged vector
	Eigen::Vector3d *toReduce, *refVector, *waiting;
	// int ids of vectors, matrix entry to transform the vector to be reduced
	int ref, red, wait, matrix_entry;
	bool reduced, face_diagonal, redirected;
	//std::cout << "Old axis (columns):\n" << axis_matrix << std::endl;

	final_transf_matrix.setIdentity();
	std::sort(axis_vector.begin(), axis_vector.end(), EigenNormCompare);
	// Start to reduce the longest with respect to the smallest
	reducing_room.push(std::make_tuple(axis_vector[2], axis_vector[1], axis_vector[0], 1));
	//std::cout << "\t|--> Reduce X with respect to (->) Y" << std::endl;
	while (!reducing_room.empty()) {
		//	Get the top element reference
		std::tuple<Eigen::Vector3d*, Eigen::Vector3d*, Eigen::Vector3d*, int> current_tuple = reducing_room.front();
		//	Remove the element
		reducing_room.pop();
		matrix_entry = 0;
		// Get pointer to vector of the top element
		toReduce = std::get<0>(current_tuple);
		refVector = std::get<1>(current_tuple);
		waiting = std::get<2>(current_tuple);

		ref = mapPointerToIndex[refVector];
		red = mapPointerToIndex[toReduce];
		wait = mapPointerToIndex[waiting];

		//Redirect the axis to converge to type 1 or type 2 cell (all acute or all obtuse angles)	
		transf_matrix.setIdentity();
		redirected = redirectAxis(axis_vector, mapPointerToIndex, transf_matrix);
		final_transf_matrix = final_transf_matrix * transf_matrix;

		reduced = reduceVector_Niggli(toReduce, refVector, waiting, face_diagonal, &matrix_entry);
		//std::cout << face_diagonal << std::endl;
		//reduced = reduceVector(toReduce, refVector, &shortest, &matrix_entry);

		// Update the transformation matrix
		if (reduced) {
			transf_matrix.setIdentity();
			// Negate the entry to transform the motif
			transf_matrix(red, ref) = -matrix_entry;
			final_transf_matrix = final_transf_matrix * transf_matrix;
		}

		// Save information about cell that is being reduced or not
		total_reduced = total_reduced || reduced || redirected;
#ifdef DEBUG
		if (reduced) { PRINTREDUCEINFO("Check Y", red, ref, -matrix_entry); }
		else { PRINTREDUCEINFO("Check N", red, ref, ""); }
#endif

		std::sort(axis_vector.begin(), axis_vector.end(), EigenNormCompare);
		//	SWITCH: Check the pointers of the first and the second and check the integer code:
		//		0: No push -> just check the second condition of a couple
		//		1: check 1st couple of vectors
		//		2: check 2nd couple of vectors
		//		3: check 3rd couple of vectors
		//		4: if fractional_proj == +-0.5 and came from the previous 0 case -> do not push
		switch (std::get<3>(current_tuple)) {
		case 0: {
			// No push -> just check the second condition of a couple
			// face_diagonal avoid pushing a case where reduction loop never stops
			// The opposite couple have been reduced and it was in face diagonal case -> Do not push
			if (reduced) {
				// A vector has been reduced, so remove the next couple -> it is invalid
				if (!reducing_room.empty()) reducing_room.pop();
				reducing_room.push(std::make_tuple(axis_vector[2], axis_vector[1], axis_vector[0], 1));
			}
			else {
				// Couple fixed -> no push
			}
			break;
		}
		case 1: {
			if (reduced) {
				//PRINTSTEPINFO(1, ref, refVector, red, toReduce, wait, waiting);
				reducing_room.push(std::make_tuple(axis_vector[2], axis_vector[1], axis_vector[0], 1));
			}
			else {
				// First couple -> Check the viceversa and push the shortest one with a new one
				reducing_room.push(std::make_tuple(axis_vector[1], axis_vector[2], axis_vector[0], 0));
				reducing_room.push(std::make_tuple(axis_vector[2], axis_vector[0], axis_vector[1], 2));
			}
			break;
		}
		case 2: {
			if (reduced) {
				//PRINTSTEPINFO(2, ref, refVector, red, toReduce, wait, waiting);
				reducing_room.push(std::make_tuple(axis_vector[2], axis_vector[1], axis_vector[0], 1));
			}
			else {
				// Second couple -> Check the viceversa and push the shortest one with the last one
				//PRINTSTEPINFO(2, red, toReduce, wait, waiting, ref, refVector);
				reducing_room.push(std::make_tuple(axis_vector[0], axis_vector[2], axis_vector[1], 0));
				reducing_room.push(std::make_tuple(axis_vector[1], axis_vector[0], axis_vector[2], 3));
			}
			break;
		}
		case 3: {
			// Last step: first 2 couple have been fixed -> check the last combination
			if (reduced) {
				//PRINTSTEPINFO(3, ref, refVector, red, toReduce, wait, waiting);
				// Enqueue to start from the beginning to check if they need to be reduced each other
				reducing_room.push(std::make_tuple(axis_vector[2], axis_vector[1], axis_vector[0], 1));
			}
			else {
				//std::cout << "\t|--> End" << std::endl;
				// Push just to check the viceversa
				reducing_room.push(std::make_tuple(axis_vector[0], axis_vector[1], axis_vector[2], 0));
			}

			break;
		}
		default: {
			break;
		}
		}

		// Check the body-diagonal only after the reduction of every vector
		if (reducing_room.empty()) {
			Eigen::Vector3d body_diagonal = Eigen::Vector3d(*(axis_vector[0]) + *(axis_vector[1]) + *(axis_vector[2]));
			//std::sort(axis_vector.begin(), axis_vector.end(), EigenNormCompare);
			Eigen::Vector3d b_face_diagonal = *(axis_vector[0]) + *(axis_vector[2]),
				a_face_diagonal = *(axis_vector[1]) + *(axis_vector[2]);
			// Guard 1: || a + b + c || < || c || or
			// The body-diagonal is smaller then the longest vector -> replace c with the body-diagonal
			// Guard 2: || a + b + c || == || c || and || a + c || > || b + c || 
			// body diagonal is equal to longest vector -> check the face-diagonals with c
			if ((body_diagonal.norm() < axis_vector[2]->norm()) ||
				((body_diagonal.norm() == axis_vector[2]->norm()) && (b_face_diagonal.norm() > a_face_diagonal.norm()))) {

				// Replace c with the body-diagonal
				std::map<Eigen::Vector3d*, int>::iterator pointer_int_it = mapPointerToIndex.find(axis_vector[2]);
				int axis_id = pointer_int_it->second;

				(*axis_vector[2])[0] = body_diagonal[0];
				(*axis_vector[2])[1] = body_diagonal[1];
				(*axis_vector[2])[2] = body_diagonal[2];

				// The 3rd vector is the 3rd + 2nd + 1st (row 2 = 1 1 1) for the axis
				// Points are transformed by row 2 = (-1,-1,1)
				transf_matrix.setIdentity();
				Eigen::Vector3d body_diag_vector = Eigen::Vector3d(-1.f, -1.f, -1.f);
				body_diag_vector[axis_id] = 1.f;
				transf_matrix.row(axis_id) = body_diag_vector;
				final_transf_matrix = final_transf_matrix * transf_matrix;
				std::sort(axis_vector.begin(), axis_vector.end(), EigenNormCompare);
				reducing_room.push(std::make_tuple(axis_vector[2], axis_vector[1], axis_vector[0], 1));
				// If the diagonal changes it stands for "reduced"
				total_reduced = true;
#ifdef DEBUG
				PRINTREDUCEINFO("BD", axis_id, "Body-Diagonal", "");
#endif
			}
		}
	}

	// Update cell parameters and the transformation matrix
	for (auto v : mapPointerToIndex) {
		new_axis.col(v.second) = *(v.first);
	}
	cell_parameters[0] = new_axis.col(0).norm();
	cell_parameters[1] = new_axis.col(1).norm();
	cell_parameters[2] = new_axis.col(2).norm();

	cell_parameters[3] = acos(new_axis.col(1).dot(new_axis.col(2)) / (new_axis.col(1).norm() * new_axis.col(2).norm())) * 180.f / M_PI;
	cell_parameters[4] = acos(new_axis.col(0).dot(new_axis.col(2)) / (new_axis.col(0).norm() * new_axis.col(2).norm())) * 180.f / M_PI;
	cell_parameters[5] = acos(new_axis.col(0).dot(new_axis.col(1)) / (new_axis.col(0).norm() * new_axis.col(1).norm())) * 180.f / M_PI;

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			transform(i, j) = final_transf_matrix(i, j);
		}
	}

	return new_axis;
	//std::cout << "New axis (columns):\n" << new_axis << std::endl;
	//std::cout << "Transformation matrix:\n" << transform << std::endl;
}
