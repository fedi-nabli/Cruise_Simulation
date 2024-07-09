/////////////////////////////////////////////////////////////
// FILE: 'Matrix.cpp'
// 'Matrix' class member functions
// 
// Module utility functions:
//	mat2tr
//	mat3tr
//	cadtge
//	cadtei
//	cadsph
//	cadtbv
//	cadine
//	sine
//	angle
// 
// Table look-up, classes 'Table' and 'Datadeck'
// Integrations
// US76 Atmosphere
// 
// Created on 07/07/2024 by Fedi Nabli
/////////////////////////////////////////////////////////////

#include "Includes/Utility.hpp"

#include <cmath>
#include <iostream>

/////////////////////////////////////////////////////////
/////////////// 'Matrix' member functions ///////////////
/////////////////////////////////////////////////////////

// Constructors
Matrix::Matrix(){}

Matrix::Matrix(int row_size, int col_size)
{
	m_NumRow = row_size;
	m_NumCol = col_size;
	m_Pbody = NULL;

	// Allocating memory
	m_NumElem = row_size * col_size;
	try
	{
		m_Pbody = new double[m_NumElem];
	}
	catch (std::bad_alloc xa)
	{
		std::cerr << "*** Error: Matrix memory allocation failed ***\n";
		std::exit(1);
	}

	// Initializing array to zero
	for (int i = 0; i < m_NumElem; i++)
		*(m_Pbody + i) = 0;
}

Matrix::Matrix(const Matrix& MAT)
{
	m_NumRow = MAT.m_NumRow;
	m_NumCol = MAT.m_NumCol;
	m_NumElem = MAT.m_NumElem;

	try
	{
		m_Pbody = new double[m_NumElem];
	}
	catch (std::bad_alloc xa)
	{
		std::cerr << "*** Error: Matrix memory allocation failed ***\n";
		std::exit(1);
	}

	// Copying
	for (int i = 0; i < m_NumElem; i++)
		*(m_Pbody + i) = (*(MAT.m_Pbody + i));
}

// Deconstructors
Matrix::~Matrix()
{
	delete[] m_Pbody;
}

// Printing Matrix to console
void Matrix::print()
{
	double* pMem = m_Pbody;

	// Outside loop rows, inside loop columns
	for (int i = 0; i < m_NumRow; i++)
	{
		for (int j = 0; j < m_NumCol; j++)
		{
			std::cout << *m_Pbody << "\t";
			m_Pbody++;
		}
		std::cout << '\n';
	}

	// Resetting pointer
	m_Pbody = pMem;
	std::cout << "\n\n";
}

// Absolute value of vector
double Matrix::absolute()
{
	if (m_NumRow > 1 && m_NumCol > 1)
	{
		std::cerr << "*** Warning: not a vector 'Matrix::absolute()' ***";
	}

	double ret = 0.0;
	for (int i = 0; i < m_NumElem; i++)
		ret += (*(m_Pbody + i)) * (*(m_Pbody + i));
	ret = std::sqrt(ret);

	return ret;
}

// Adjoint matrix (same as determinant procedure however the matrix element
// is NOT multiplied into each cofactor)
Matrix Matrix::adjoint()
{
	if (!(m_NumRow == m_NumCol))
	{
		std::cerr << "*** Error: matrix not square 'Matrix::adjoint()' ***\n";
		std::exit(1);
	}

	if ((m_NumRow == 1) && (m_NumCol == 1))
	{
		std::cerr << "*** Error: only one element 'Matrix::adjoint()' ***\n";
		std::exit(1);
	}

	Matrix RESULT(m_NumCol, m_NumCol);

	for (int i = 0; i < m_NumElem; i++)
	{
		// row #
		int row = i / m_NumCol + 1;
		// column #
		int col = i % m_NumCol + 1;

		if (((row + col) % 2) == 0)
			*(RESULT.m_Pbody + i) = sub_matrix(row, col).determinant();
		else
			*(RESULT.m_Pbody + i) = (-1.0) * sub_matrix(row, col).determinant();
	}

	return RESULT.trans();
}

// Assigns a value to a matrix element (offset)
void Matrix::assign_loc(const int& r, const int& c, const double& val)
{
	if (r > m_NumRow - 1 || c > m_NumCol - 1)
	{
		std::cerr << "*** Error: location outside array 'Matrix::assign_loc()' ***\n";
		std::exit(1);
	}

	// Assign a value
	int offset = m_NumCol * (r) + c;
	*(m_Pbody + offset) = val;
}

// Builds a 3x1 vector from parameters
Matrix& Matrix::build_vec3(const double& v1, const double& v2, const double& v3)
{
	m_NumRow = 3;
	m_NumCol = 1;
	*m_Pbody = v1;
	*(m_Pbody + 1) = v2;
	*(m_Pbody + 2) = v3;

	return *this;
}

Matrix& Matrix::build_mat33(const double& v11, const double& v12, const double& v13,
	const double& v21, const double& v22, const double& v23,
	const double& v31, const double& v32, const double& v33)
{
	m_NumRow = 3;
	m_NumCol = 3;
	*m_Pbody = v11;
	*(m_Pbody + 1) = v12;
	*(m_Pbody + 2) = v13;
	*(m_Pbody + 3) = v21;
	*(m_Pbody + 4) = v22;
	*(m_Pbody + 5) = v23;
	*(m_Pbody + 6) = v31;
	*(m_Pbody + 7) = v32;
	*(m_Pbody + 8) = v33;

	return *this;
}

// Returns vector of row #
Matrix Matrix::row_vec(const int& row)
{
	if (row <= 0 || row > m_NumRow)
	{
		std::cerr << "*** Error: row outside array 'Matrix::row_vec()' ***\n";
		std::exit(1);
	}

	Matrix RESULT(m_NumCol, 1);

	for (int i = 0; i < m_NumCol; i++)
	{
		int offset = (row - 1) * m_NumCol + i;
		*(RESULT.m_Pbody + i) = (*(m_Pbody + offset));
	}

	return RESULT;
}

// Returns column vector of column #
Matrix Matrix::col_vec(const int& col)
{
	if (col <= 0 || col > m_NumCol)
	{
		std::cerr << "*** Error: column outside array 'Matrix::col_vec()' ***\n";
		std::exit(1);
	}

	Matrix RESULT(m_NumRow, 1);

	for (int i = 0; i < m_NumRow; i++)
	{
		int offset = i * m_NumCol + col - 1;
		*(RESULT.m_Pbody + i) = (*(m_Pbody + offset));
	}

	return RESULT;
}

double Matrix::determinant()
{
	if (!(m_NumRow == m_NumCol))
	{
		std::cerr << "*** Error: matrix not square 'Matrix::determinant()' ***\n";
		std::exit(1);
	}

	double result = 0.0;

	// Base case of a single matrix element
	if ((m_NumCol == 1) && (m_NumRow == 1))
		return *m_Pbody;

	else if ((m_NumCol == 2) && (m_NumRow == 2))
		return (*m_Pbody) * (*(m_Pbody + 3)) - (*(m_Pbody + 1)) * (*(m_Pbody + 2));

	else
	{
		for (int j = 0; j < m_NumCol; j++)
		{
			// Use cofactors and submatricies to finish for NxN
			if ((j % 2) == 0)
			{
				// Odd column (numbered)
				result += sub_matrix(1, j + 1).determinant() * (*(m_Pbody + j));
			}
			else
			{
				// Even column
				result += (-1.0) * sub_matrix(1, j + 1).determinant() * (*(m_Pbody + j));
			}
		}
	}

	return result;
}

// Returns NxN diagonal matrix from Nx1 vector
Matrix Matrix::diamat_vec()
{
	if (m_NumCol != 1)
	{
		std::cerr << "*** Error: not a vector 'Matrix::diamet_vec()' ***\n";
		std::exit(1);
	}

	Matrix RESULT(m_NumRow, m_NumCol);
	for (int i = 0; i < m_NumRow; i++)
	{
		int offset = i * m_NumRow + i;
		*(RESULT.m_Pbody + offset) = (*(m_Pbody + i));
	}

	return RESULT;
}

// Returns Nx1 diagonal vector from NxN matrix
Matrix Matrix::diavec_mat()
{
	if (!(m_NumRow == m_NumCol))
	{
		std::cerr << "*** Error: matrix not square 'Matrix::diavec_mat()' ***\n";
		std::exit(1);
	}

	Matrix RESULT(m_NumRow, 1);
	for (int i = 0; i < m_NumRow; i++)
	{
		int offset = i * m_NumRow + i;
		*(RESULT.m_Pbody + i) = (*(m_Pbody + offset));
	}

	return RESULT;
}

// Dimensions a matrix of size row x col
// Only used to initliaze arrays in class 'Variable'
void Matrix::dimension(int row, int col)
{
	m_NumRow = row;
	m_NumCol = col;
	m_Pbody = NULL;

	// Allocating memory
	m_NumElem = row * col;
	try
	{
		m_Pbody = new double[m_NumElem];
	}
	catch (std::bad_alloc xa)
	{
		std::cerr << "*** Error: memory allocation failed 'Matrix::dimension()' ***\n";
		std::exit(1);
	}

	// Initilaizing array to zero
	for (int i = 0; i < m_NumElem; i++)
		*(m_Pbody + i) = 0.;
}

// Returns the number of rows of matrix MAT
int Matrix::get_rows()
{
	return m_NumRow;
}

// Returns the number of columns of matrix MAT
int Matrix::get_cols()
{
	return m_NumCol;
}

// Returns offset-index fiven row# and col#
int Matrix::get_index(const int& row, const int& col)
{
	int index;
	index = (row - 1) * m_NumCol + col - 1;
	return index;
}

// Returns the value at offset-row 'r' and offset-col 'c' of MAT
double Matrix::get_loc(const int& r, const int& c)
{
	if ((r < m_NumRow) && (c < m_NumCol))
		return *(m_Pbody + r * m_NumCol + c);
	else
	{
		std::cerr << "*** Error: invalid matrix location 'Matrix::get_loc()' ***\n";
		std::exit(1);
	}
}

double* Matrix::get_pbody()
{
	return m_Pbody;
}

// Calculates Cartesian vector from polar coordinates
// |V1|				  | cos(elevation) * cos(azimuth) |
// |V2| = magnitude * | cos(elevation) * sin(azimuth) |
// |V3|				  |		  -sin(elevation)         |
Matrix& Matrix::cart_from_pol(const double& magnitude, const double& azimuth, const double& elevation)
{
	*m_Pbody = magnitude * (std::cos(elevation) * std::cos(azimuth));
	*(m_Pbody + 1) = magnitude * (std::cos(elevation) * std::sin(azimuth));
	*(m_Pbody + 2) = magnitude * (std::sin(elevation) * (-1.0));

	return *this;
}

// Returns polar from cartesian coordinates
// magnitude = POLAR(0,0) = |V|
// azimuth	 = POLAR(1,0) = atan2(V2, V1)
// elevation = POLAR(2,0) = atan2(-V3, sqrt(V1^2+V2^2))
Matrix Matrix::pol_from_cart()
{
	double d = 0.0;
	double azimuth = 0.0;
	double elevation = 0.0;
	double denom;

	Matrix POLAR(3, 1);

	double v1 = (*m_Pbody);
	double v2 = (*(m_Pbody + 1));
	double v3 = (*(m_Pbody + 2));

	d = std::sqrt(v1 * v1 + v2 * v2 + v3 * v3);
	azimuth = std::atan2(v2, v1);

	denom = std::sqrt(v1 * v1 + v2 * v2);
	if (denom > 0.)
		elevation = std::atan2(-v3, denom);
	else
	{
		if (v3 > 0)
			elevation = -PI / 2.;
		if (v3 < 0)
			elevation = PI / 2.;
		if (v3 == 0)
			elevation = 0.;
	}

	*POLAR.m_Pbody = d;
	*(POLAR.m_Pbody + 1) = azimuth;
	*(POLAR.m_Pbody + 2) = elevation;

	return POLAR;
}

///////////////////////////////////////////////////////////////////////////////
// Bi-variate ellipse
// calculating major and minor semi-axes of ellipse and rotation angle 
//    from the symmetrical pos semi-definite MAT(2x2) matrix
// coordinate axes orientation:
//          ^ 1-axis
//          |
//          |
//          |---> 2-axis
//
// angle is measured from 1st coordinate axis to the right
//
// major_semi_axis = ELLIPSE.get_loc(0,0);
// minor_semi_axis = ELLIPSE.get_loc(1,0);
// angle      = ELLIPSE.get_loc(2,0);
///////////////////////////////////////////////////////////////////////////////
Matrix Matrix::ellipse()
{
	Matrix ELLIPSE(3, 1);

	double dum = 0;
	double dum1 = 0;
	double dum2 = 0;
	double ama = 0;
	double ami = 0;
	double phi = 0;
	double ak1 = 0;
	double ak2 = 0;

	Matrix X1V(2, 1); // Major principal axes of ellipse
	Matrix X2V(2, 1); // Minor principal axes of ellipse

	double a11 = *m_Pbody;
	double a12 = *(m_Pbody + 1);
	double a22 = *(m_Pbody + 3);
	double a1122 = a11 + a22;
	double aq1122 = a1122 * a1122;
	dum1 = aq1122 - 4. * (a11 * a22 - a12 * a12);
	if (dum1 >= 0.)
		dum2 = std::sqrt(dum1);

	// Major and minor semi-axes of ellipse
	ama = (a1122 + dum2) / 2.;
	ami = (a1122 - dum2) / 2.;
	ELLIPSE.assign_loc(0, 0, ama);
	ELLIPSE.assign_loc(1, 0, ami);
	if (ama == ami)
		return ELLIPSE;

	// Angle of orientation of major axis wrt first peincipal axis
	if (a11 - ama != 0.)
	{
		dum1 = -a12 / (a11 - ama);
		ak1 = std::sqrt(1. / (1. + dum1 * dum1));
		X1V.assign_loc(0, 0, dum1 * ak1);
		X1V.assign_loc(1, 0, ak1);
		dum = dum1 * ak1;
		if (std::fabs(dum) > 1.)
			dum = 1. * sign(dum);
		phi = std::acos(dum);
		ELLIPSE.assign_loc(2, 0, phi);
	}
	else
	{
		dum1 = -a12 / (a22 - ama);
		ak1 = std::sqrt(1. / (1. + dum1 * dum1));
		X1V.assign_loc(0, 0, ak1);
		X1V.assign_loc(1, 0, dum1 * dum1);
		if (std::fabs(ak1) > 1.)
			ak1 = 1. * sign(ak1);
		phi = std::acos(ak1);
		ELLIPSE.assign_loc(2, 0, phi);
	}

	// Second principal axis - not used
	if (a11 - ami != 0.)
	{
		dum2 = -a12 / (a11 - ami);
		ak2 = std::sqrt(1. / (1. + dum2 * dum2));
		X2V.assign_loc(0, 0, dum2 * ak2);
		X2V.assign_loc(1, 0, ak2);
	}
	else
	{
		dum2 = -a12 / (a22 - ami);
		ak2 = std::sqrt(1. / (1. + dum2 * dum2));
		X2V.assign_loc(0, 0, ak2);
		X2V.assign_loc(1, 0, ak2);
	}

	return ELLIPSE;
}

// Builds a square identity matrix of object 'Matrix MAT'
Matrix& Matrix::identity()
{
	if (m_NumRow == m_NumCol)
	{
		for (int r = 0; r < m_NumRow; r++)
			*(m_Pbody + r * m_NumRow + r) = 1.;
	}
	else
	{
		std::cerr << "*** Error: matrix not square 'Matrix::identity()' ***\n";
		std::exit(1);
	}

	return *this;
}

// Returns the inverse of a square matrix MAT
// Inversion: Inverse = (1/det(A)) * Adj(A)
Matrix Matrix::inverse()
{
	if (m_NumCol != m_NumRow)
	{
		std::cerr << "*** Error: not a square matrix 'Matrix::inverse()' ***\n";
		std::exit(1);
;	}

	Matrix RESULT(m_NumRow, m_NumCol);
	double d = 0.;

	d = determinant();
	if (d == 0.)
	{
		std::cerr << "*** Error: singular! 'Matrix::inverse()' ***\n";
		std::exit(1);
	}

	d = 1. / d;
	RESULT = adjoint();
	RESULT = RESULT * d;

	return RESULT;
}

// Returns 3x3 marix row-size from 9x1 vector
Matrix Matrix::mat33_vec9()
{
	if (!(m_NumRow == 9 && m_NumCol == 1))
	{
		std::cerr << "*** Error: vector not 9 x 1 'Matrix::mat33_vec9()' ***\n";
		std::exit(1);
	}

	Matrix RESULT(3, 3);
	for (int i = 0; i < 9; i++)
	{
		*(RESULT.m_Pbody + i) = *(m_Pbody + i);
	}

	return RESULT;
}

// Forms matrix MAT with all elementes '1.'
Matrix& Matrix::ones()
{
	for (int r = 0; r < m_NumElem; r++)
		*(m_Pbody + r) = 1.;

	return *this;
}

// Returns the skew-symmetric matrix MAT from a 3-dim vector VEC
//		| 0 -c  b|	   |a|
//		| c  0 -a| <-- |b|
//		|-b  a  0|	   |c|
Matrix Matrix::skew_sym()
{
	Matrix RESULT(3, 3);
	
	// Check for proper dimensions
	if (m_NumCol != 1 || m_NumRow != 3)
	{
		std::cerr << "*** Error: not a 3x1 column vector 'Matrix::skew_sym()' ***\n";
		std::exit(1);
	}

	*(RESULT.m_Pbody + 5) = -(*m_Pbody);
	*(RESULT.m_Pbody + 7) = (*m_Pbody);
	*(RESULT.m_Pbody + 2) = (*(m_Pbody + 1));
	*(RESULT.m_Pbody + 6) = -(*(m_Pbody + 1));
	*(RESULT.m_Pbody + 1) = -(*(m_Pbody + 2));
	*(RESULT.m_Pbody + 3) = (*(m_Pbody + 2));

	return RESULT;
}

// Returns the sub matrix after 'row' and 'col' have been ommitted
Matrix Matrix::sub_matrix(const int& row, const int& col)
{
	if ((row > m_NumRow) || (col > m_NumCol))
	{
		std::cerr << "*** Error: row or column outside array 'Matrix::sub_matrix()' ***\n";
		std::exit(1);
	}

	if (row == 0 || col == 0)
	{
		std::cerr << "*** Error: row/col are numbered not offset 'Matrix::sub_matrix()' ***\n";
		std::exit(1);
	}

	// Create return matrix
	Matrix RESULT(m_NumRow - 1, m_NumCol - 1);
	
	// Start and Stop of skipping matrix elements
	int skip_start = (row - 1) * m_NumCol;
	int skip_end = skip_start + m_NumCol;

	// Intiailize  RESULT offset j
	int j = 0;

	for (int i = 0; i < m_NumElem; i++)
	{
		// Skip elements of row to be romoved
		if ((i < skip_start) || (i >= skip_end))
		{
			// Offset of column element to be removed
			int offset_col = (col - 1) + (i / m_NumCol) * m_NumCol;
			// Skip elements of col to be removed
			if (i != offset_col)
			{
				*(RESULT.m_Pbody + j) = *(m_Pbody + i);
				j++;
			}
		}
	}

	return RESULT;
}

// Returns the transpose of a matrix Aji < Aij
// (same as 'Matrix operator~')
Matrix Matrix::trans()
{
	Matrix RESULT(m_NumCol, m_NumRow);

	int i = 0; // Offset for original matrix
	int j = 0; // Offset for transposed matrix

	for (int r = 0; r < m_NumRow; r++)
	{
		for (int c = 0; c < m_NumCol; c++)
		{
			// Offset for transposed
			j = c * m_NumRow + r;
			*(RESULT.m_Pbody + j) = (*(m_Pbody + i));
			i++;
			j++;
		}
	}

	return RESULT;
}

// Returns unit vector from 3x1 vector
Matrix Matrix::univec3()
{
	Matrix RESULT(3, 1);

	// Check for proper dimensions
	if (m_NumCol != 1 || m_NumRow != 3)
	{
		std::cerr << "*** Error: not a 3x1 column vector in 'Matrix::univec3()' ***\n";
		std::exit(1);
	}

	double v1 = (*m_Pbody);
	double v2 = (*(m_Pbody + 1));
	double v3 = (*(m_Pbody + 2));
	
	double d = std::sqrt(v1 * v1 + v2 * v2 + v3 * v3);

	// if VEC is zero than the unit vector is also a zero vector
	if (d == 0)
	{
		*(RESULT.m_Pbody) = 0;
		*(RESULT.m_Pbody + 1) = 0;
		*(RESULT.m_Pbody + 2) = 0;
	}
	else
	{
		*(RESULT.m_Pbody) = v1 / d;
		*(RESULT.m_Pbody + 1) = v2 / d;
		*(RESULT.m_Pbody + 2) = v3 / d;
	}

	return RESULT;
}

// Returns 9x1 vector from 3x3 matrix
Matrix Matrix::vec9_mat33()
{
	if (!(m_NumRow == 3 && m_NumCol == 3))
	{
		std::cerr << "*** Error: matrix not 3x3 'Matrix::vec9_mat33()' ***\n";
		std::exit(1);
	}

	Matrix RESULT(9, 1);

	for (int i = 0; i < 9; i++)
		*(RESULT.m_Pbody + i) = *(m_Pbody + i);

	return RESULT;
}

// Forms a zero matrix MAT from object MAT(num_row, num_col)
Matrix& Matrix::zero()
{
	for (int i = 0; i < m_NumElem; i++)
		*(m_Pbody + i) = 0.0;

	return *this;
}

// Inequality relational operator, returns true of false
// returns true if elements differe by more than EPS
bool Matrix::operator!=(const Matrix& B)
{
	// Check dimensions
	if (m_NumCol != B.m_NumCol)
		return true;

	else if (m_NumRow != B.m_NumRow)
		return true;

	for (int i = 0; i < m_NumElem; i++)
	{
		// Check to see if values differ by more or less than EPS
		if ((*(m_Pbody + i) - (*(B.m_Pbody + i))) > EPS)
			return true;
		else if ((*(m_Pbody + i) - (*(B.m_Pbody + i))) < (-1. * EPS))
			return true;
	}

	return false;
}

// Equality relational operator
bool Matrix::operator==(const Matrix& B)
{
	// Check dimensions
	if (m_NumCol != B.m_NumCol)
		return false;

	else if (m_NumRow != B.m_NumRow)
		return false;

	for (int i = 0; i < m_NumElem; i++)
	{
		// Check to see if values differ by more or less than EPS
		if ((*(m_Pbody + i) - (*(B.m_Pbody + i))) > EPS)
			return false;
		else if ((*(m_Pbody + i) - (*(B.m_Pbody + i))) < (-1. * EPS))
			return false;
	}

	return true;
}

// Scalar multiplication operator (scalar element by element multiplicaton)
// Note: scalar must be the second operand
Matrix Matrix::operator*(const double& b)
{
	Matrix RESULT(m_NumRow, m_NumCol);

	for (int i = 0; i < m_NumElem; i++)
		*(RESULT.m_Pbody + i) = *(m_Pbody + i) * b;

	return RESULT;
}

// Multiplication operator, returns matrix product
Matrix Matrix::operator*(const Matrix& B)
{
	// create resultant matrix
	Matrix RESULT(m_NumRow, B.m_NumCol);
	int r = 0;
	int c = 0;

	// Check for proper dimensions
	if (m_NumCol != B.m_NumRow)
	{
		std::cerr << "*** Error: incompatible dimensions 'Matrix::operator*()' ***\n";
		std::exit(1);
	}

	for (int i = 0; i < RESULT.m_NumElem; i++)
	{
		r = i / B.m_NumCol;
		c = i % B.m_NumCol;
		for (int k = 0; k < m_NumCol; k++)
		{
			*(RESULT.m_Pbody + i) += *(m_Pbody + k + m_NumCol * r) * (*(B.m_Pbody + k * B.m_NumCol + c));
		}
	}

	return RESULT;
}

// Scalar multiplication assignment operator
Matrix& Matrix::operator*=(const double& b)
{
	for (int i = 0; i < m_NumElem; i++)
	{
		*(m_Pbody + i) = *(m_Pbody + i) * b;
	}

	return *this;
}

// Multiplication assignment operator
// Matrix B must be square
Matrix& Matrix::operator*=(const Matrix& B)
{
	// Create resultant matrix
	Matrix RESULT(m_NumRow, B.m_NumCol);

	// Check for proper dimensions
	if (m_NumCol != B.m_NumRow)
	{
		std::cerr << "*** Error: incompatible dimensions 'Matrix::operatror*=()' ***\n";
		std::exit(1);
	}

	// Check for squareness of B
	if (B.m_NumRow != B.m_NumCol)
	{
		std::cerr << "*** Error: second matrix is not square 'Matrix::operator*=()' ***\n";
		std::exit(1);
	}

	int i;
	for (i = 0; i < m_NumElem; i++)
	{
		int r = i / B.m_NumCol;
		int c = i % B.m_NumCol;
		for (int k = 0; k < m_NumCol; k++)
		{
			*(RESULT.m_Pbody + i) += *(m_Pbody + k + m_NumCol * r) * (*(B.m_Pbody + k + B.m_NumCol + c));
		}
	}

	m_NumCol = RESULT.m_NumCol;
	m_NumRow = RESULT.m_NumRow;
	m_NumElem = m_NumRow * m_NumCol;
	for (i = 0; i < m_NumElem; i++)
		*(m_Pbody + i) = *(RESULT.m_Pbody + i);

	return *this;
}

// Scalar addition operator
// Note scarlat must be the second operand
Matrix Matrix::operator+(const double& b)
{
	Matrix RESULT(m_NumRow, m_NumCol);

	for (int i = 0; i < m_NumElem; i++)
		*(RESULT.m_Pbody + i) = *(m_Pbody + i) + b;

	return RESULT;
}

// Addition operator, returns matrix addition
// Operands must be conformal
Matrix Matrix::operator+(const Matrix& B)
{
	Matrix RESULT(m_NumRow, m_NumCol);

	if ((m_NumCol != B.m_NumCol) || (m_NumRow != B.m_NumRow))
	{
		std::cerr << "*** Error: matrices have different dimensions 'Matrix::operator+()' ***\n";
		std::exit(1);
	}

	for (int i = 0; i < m_NumElem; i++)
		*(RESULT.m_Pbody + i) = *(m_Pbody + i) + (*(B.m_Pbody + i));

	return RESULT;
}

// Scalar addition assignment opeeator
Matrix& Matrix::operator+=(const double& b)
{
	for (int i = 0; i < m_NumElem; i++)
		*(m_Pbody + i) = *(m_Pbody + i) + b;

	return *this;
}

// Matrix addition assignment operator
Matrix& Matrix::operator+=(const Matrix& B)
{
	if ((m_NumCol != B.m_NumCol) || (m_NumRow != B.m_NumRow))
	{
		std::cerr << "*** Error: matrices have different diemnsions 'Matrix::operator+=()' ***\n";
		std::exit(1);
	}

	for (int i = 0; i < m_NumElem; i++)
		*(m_Pbody + i) = *(m_Pbody + i) + (*(B.m_Pbody + i));

	return *this;
}

// Scalar substruction operator
// Note: scalar must be the second operand
Matrix Matrix::operator-(const double& b)
{
	Matrix RESULT(m_NumRow, m_NumCol);
	
	for (int i = 0; i < m_NumElem; i++)
		*(RESULT.m_Pbody + i) = *(m_Pbody + i) - b;

	return RESULT;
}

// Substruction operator, returns new matrix
Matrix Matrix::operator-(const Matrix& B)
{
	Matrix RESULT(m_NumRow, m_NumCol);

	if ((m_NumCol != B.m_NumCol) || (m_NumRow != B.m_NumRow))
	{
		std::cerr << "*** Error: matrices have different dimensons 'Matrix::operator-()' ***\n";
		std::exit(1);
	}

	for (int i = 0; i < m_NumElem; i++)
		*(RESULT.m_Pbody + i) = *(m_Pbody + i) - (*(B.m_Pbody + i));

	return RESULT;
}

// Scalar substruction assignment operator
Matrix& Matrix::operator-=(const double& b)
{
	for (int i = 0; i < m_NumElem; i++)
		*(m_Pbody + i) = *(m_Pbody + i) - b;

	return *this;
}

// Matrix substruction assignment operator
Matrix& Matrix::operator-=(const Matrix& B)
{
	if ((m_NumCol != B.m_NumCol) || (m_NumRow != B.m_NumRow))
	{
		std::cerr << "*** Error: matrices have different dimensions 'Matrix::operator-=()' ***\n";
		std::exit(1);
	}

	for (int i = 0; i < m_NumElem; i++)
		*(m_Pbody + i) = *(m_Pbody + i) - (*(B.m_Pbody + i));

	return *this;
}

// Assignment operator (deep copy)
Matrix& Matrix::operator=(const Matrix& B)
{
	if ((m_NumCol != B.m_NumCol) || (m_NumRow != B.m_NumRow))
	{
		std::cerr << "*** Error: incompatible dimensions 'Matrix::operator=()' ***\n";
		std::exit(1);
	}

	delete[] m_Pbody;
	m_NumElem = B.m_NumElem;
	m_NumRow = B.m_NumRow;
	m_NumCol = B.m_NumCol;
	m_Pbody = new double[m_NumElem];

	for (int i = 0; i < m_NumElem; i++)
		*(m_Pbody + i) = (*(B.m_Pbody + i));

	return *this;
}

// Returns the component(i) from vector VEC[i] or assigns a value to component VEC[i]
double& Matrix::operator[](const int& r)
{
	if ((r < m_NumRow) && (m_NumCol == 1))
		return *(m_Pbody + r);
	else
	{
		std::cerr << "*** Error: invalid matrix location 'Matrix::operator[]()' ***\n";
		std::system("pause");
		std::exit(1);
	}
}

// Scalar product operator (any combination of row or column vectors)
double Matrix::operator^(const Matrix& B)
{
	// Initialize the result
	double result = 0.0;

	// Check diemensions
	bool one = false;
	bool dim = false;

	// True if both arrays have dimension '1'
	if ((m_NumRow == 1 || m_NumCol == 1) && (B.m_NumRow == 1 || B.m_NumCol == 1))
		one = true;

	// True if both arrays have at least one equal dimension
	if ((m_NumRow == B.m_NumRow || m_NumRow == B.m_NumCol) && (m_NumCol == B.m_NumCol || m_NumCol == B.m_NumRow))
		dim = true;

	if (!one || !dim)
	{
		std::cerr << "*** Error: incompatible dimensions 'Matrix::operator^()' ***\n";
		std::exit(1);
	}

	for (int i = 0; i < m_NumElem; i++)
		result += *(m_Pbody + i) * (*(B.m_Pbody + i));

	return result;
}

// Alternate transpose Aji < Aij
Matrix Matrix::operator~()
{
	Matrix RESULT(m_NumCol, m_NumRow);

	int i = 0; // Offset for original matrix
	int j = 0; // Offset for transposed matrix

	for (int r = 0; r < m_NumRow; r++)
	{
		for (int c = 0; c < m_NumCol; c++)
		{
			// Offset for transposed
			j = c * m_NumRow + r;
			*(RESULT.m_Pbody + j) = *(m_Pbody + i);
			i++;
			j++;
		}
	}

	return RESULT;
}


//////////////////////////////////////////////////////////
//////////////// Module utility functions ////////////////
//////////////////////////////////////////////////////////

// Returns the T.M. of psivg -> thtvg sequence
Matrix mat2tr(const double& psivg, const double& thtvg)
{
	Matrix AMAT(3, 3);

	AMAT.assign_loc(0, 2, -std::sin(thtvg));
	AMAT.assign_loc(1, 0, -std::sin(psivg));
	AMAT.assign_loc(1, 1, std::cos(psivg));
	AMAT.assign_loc(2, 2, std::cos(thtvg));
	AMAT.assign_loc(0, 0, (AMAT.get_loc(2, 2) * AMAT.get_loc(1, 1)));
	AMAT.assign_loc(0, 1, (-AMAT.get_loc(2, 2) * AMAT.get_loc(1, 0)));
	AMAT.assign_loc(2, 0, (-AMAT.get_loc(0, 2) * AMAT.get_loc(1, 1)));
	AMAT.assign_loc(2, 1, (AMAT.get_loc(0, 2) * AMAT.get_loc(1, 0)));
	AMAT.assign_loc(1, 2, 0.0);

	return AMAT;
}

// Returns the Euler T.M. of the psi -> tht -> phi sequence
Matrix mat3tr(const double& psi, const double& tht, const double& phi)
{
	double spsi = std::sin(psi);
	double cpsi = std::cos(psi);
	double stht = std::sin(tht);
	double ctht = std::cos(tht);
	double sphi = std::sin(phi);
	double cphi = std::cos(phi);

	Matrix AMAT(3, 3);

	AMAT.assign_loc(0, 0, cpsi * ctht);
	AMAT.assign_loc(1, 0, cpsi * stht * sphi - spsi * cphi);
	AMAT.assign_loc(2, 0, cpsi * stht * cphi + spsi * sphi);
	AMAT.assign_loc(0, 1, spsi * ctht);
	AMAT.assign_loc(1, 1, spsi * stht * sphi + cpsi * cphi);
	AMAT.assign_loc(2, 1, spsi * stht * cphi - cpsi * sphi);
	AMAT.assign_loc(0, 2, -stht);
	AMAT.assign_loc(1, 2, ctht * sphi);
	AMAT.assign_loc(2, 2, ctht * cphi);

	return AMAT;
}

// Returns the T.M. of geographic wrt earth coordinates
Matrix cadtge(double lon, double lat)
{
	Matrix AMAT(3, 3);

	double clon(0);
	double slon(0);
	double clat(0);
	double slat(0);

	clon = std::cos(lon);
	slon = std::sin(lon);
	clat = std::cos(lat);
	slat = std::sin(lat);

	AMAT.assign_loc(0, 0, (-slat * clon));
	AMAT.assign_loc(0, 1, (-slat * slon));
	AMAT.assign_loc(0, 2, clat);
	AMAT.assign_loc(1, 0, -slon);
	AMAT.assign_loc(1, 1, clon);
	AMAT.assign_loc(1, 2, 0.0);
	AMAT.assign_loc(2, 0, (-clat * clon));
	AMAT.assign_loc(2, 1, (-clat * slon));
	AMAT.assign_loc(2, 2, -slat);

	return AMAT;
}

// Returns the T.M. of earth wrt inertial coordinates
Matrix cadtei(double simulation_time)
{
	double xi(0);
	double sxi(0);
	double cxi(0);

	Matrix TEI(3, 3);

	xi = WEII3 * simulation_time;
	sxi = std::sin(xi);
	cxi = std::cos(xi);

	TEI.identity();
	TEI.assign_loc(0, 0, cxi);
	TEI.assign_loc(0, 1, sxi);
	TEI.assign_loc(1, 0, -sxi);
	TEI.assign_loc(1, 1, cxi);

	return TEI;
}

// Returns lon, lat, alt from inertial displacement vector
Matrix cadsph(Matrix SBIE)
{
	double dum4(0);
	double alamda(0);
	double x(0), y(0), z(0);

	x = SBIE.get_loc(0, 0);
	y = SBIE.get_loc(1, 0);
	z = SBIE.get_loc(2, 0);

	double dbi;
	double alt;
	double lat;
	double lon;

	Matrix RESULT(3, 1);

	// Latitude
	dbi = std::sqrt(x * x + y * y + z * z);
	lat = std::asin((z) / dbi);

	// Altitude
	alt = dbi - REARTH;

	// Longitude
	dum4 = std::asin(y / std::sqrt(x * x + y * y));

	// Resolving the multi-valued arcsin function
	if ((x >= 0) && (y >= 0))
	{
		alamda = dum4; // Quadrant I
	}

	if ((x < 0) && (y >= 0))
	{
		alamda = 180 * RAD - dum4; // Quadrant II
	}

	if ((x < 0) && (y < 0))
	{
		alamda = 180 * RAD - dum4; // Quadrant III
	}

	if ((x >= 0) && (y < 0))
	{
		alamda = 360 * RAD + dum4; // Quadrant IV
	}

	lon = alamda;

	if (lon > (180 * RAD))
	{
		lon = -(360 * RAD - lon); // east positive, west negative
	}

	RESULT.assign_loc(0, 0, lon);
	RESULT.assign_loc(1, 0, lat);
	RESULT.assign_loc(2, 0, alt);

	return RESULT;
}

// Returns the T.M. of body wrt velocity coordinates (3 DoF)
// Suitable for 3 DoF back-to-turn simulations only
Matrix cadtbv(double phi, double alpha)
{
	Matrix AMAT(3, 3);

	double sphi = std::sin(phi);
	double cphi = std::cos(phi);
	double salpha = std::sin(alpha);
	double calpha = std::cos(alpha);
	
	AMAT.assign_loc(0, 0, calpha);
	AMAT.assign_loc(0, 1, sphi * salpha);
	AMAT.assign_loc(0, 2, -cphi * salpha);
	AMAT.assign_loc(1, 0, 0.0);
	AMAT.assign_loc(1, 1, cphi);
	AMAT.assign_loc(1, 2, sphi);
	AMAT.assign_loc(2, 0, salpha);
	AMAT.assign_loc(2, 1, -sphi * calpha);
	AMAT.assign_loc(2, 2, cphi * calpha);

	return AMAT;
}

// Returns inertial coordinates from longitude, latitude and altitude
Matrix cadine(double lon, double lat, double alt, double time)
{
	Matrix VEC(3, 1);

	double rad = alt + REARTH;
	double cel_lon = lon + WEII3 * time;
	double clat = std::cos(lat);
	double slat = std::sin(lat);
	double clon = std::cos(cel_lon);
	double slon = std::sin(cel_lon);

	VEC.assign_loc(0, 0, rad * clat * clon);
	VEC.assign_loc(1, 0, rad * clat * slon);
	VEC.assign_loc(2, 0, rad * slat);

	return VEC;
}

// Returns the sign of the function
int sign(double variable)
{
	int sign(0);

	if (variable < 0)
		sign = -1;
	if (variable >= 0)
		sign = 1;

	return sign;
}

// Returns the angle between two 3x1 vectors
double angle(Matrix VEC1, Matrix VEC2)
{
	double argument(0);
	double scalar = VEC1 ^ VEC2;
	double abs1 = VEC1.absolute();
	double abs2 = VEC2.absolute();

	double dum = abs1 * abs2;

	if ((abs1 * abs2) > EPS)
		argument = scalar / dum;
	else
		argument = 1.;

	if (argument > 1.)
		argument = 1.;
	if (argument < -1.)
		argument = -1.;

	return std::acos(argument);
}

int main()
{
	Matrix TEST(3, 3);
	TEST.print();

	TEST.identity();
	TEST.print();

	TEST.ones();
	TEST.print();

	Matrix RESULT = cadsph(TEST);
	RESULT.print(),

	TEST.zero();
	TEST.print();

	return 0;
}