/////////////////////////////////////////////////////////////
// FILE: 'Utility.hpp'
// Declares the classes: 'Matrix', 'Table' and 'Datadeck'.
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

#pragma once

#include "GlobalConstants.hpp"

#include <cmath>
#include <iostream>

/// <summary>
/// 'Matrix' member functions
/// One dimensional and two dimensional arrays of any dimension of type 'double'
/// Class 'Matrix'
///	 dynamically allocated matrix size
///  pointer arithmetic
/// </summary>

class Matrix
{
public:
	// Default constructor
	Matrix();
	// Overloaded constructor
	Matrix(int row_size, int col_size);
	// Copy constructor
	Matrix(const Matrix& MAT);
	// Deconstructor
	~Matrix();

	// Printing matrix array to console
	void print();

	// Returns absolute value of a vector VEC
	double absolute();

	// Returns the adjoint
	Matrix adjoint();

	// Assigns a value tp a matrix element (offset)
	void assign_loc(const int& r, const int& c, const double& val);

	// Builds a 3x1 vector from three parameters
	Matrix& build_vec3(const double& v1, const double& v2, const double& v3);

	Matrix& build_mat33(const double& v11, const double& v12, const double& v13,
		const double& V21, const double& v22, const double& v33,
		const double& v31, const double& v32, const double& V33);

	// Returns row vector of row #
	Matrix row_vec(const int& row);

	// Returns column vector of column #
	Matrix col_vec(const int& col);

	// Returns the determinant
	double determinant();

	// Returns NxN diagonal matrix from Nx1 vector
	Matrix diamat_vec();

	// Returns Nx1 diagonal vector from NxN matrix
	Matrix diavec_mat();

	// Dimensions a matrix of size row x col
	void dimension(int row, int col);

	// Returns the number of rows in the matrix
	int get_rows();

	// Returns the number of columns of matrix MAT
	int get_cols();

	// Returns offset index given row# and col#
	int get_index(const int& row, const int& col);

	// Returns the value at offset-row 'r' offset-col 'c' of MAT
	double get_loc(const int& r, const int& c);

	// Returnds the pointer to MAT
	double* get_pbody();

	// Calculates Cartesian vector from polar coordinates
	// |V1|				  | cos(elevation) * cos(azimuth) |
	// |V2| = magnitude * | cos(elevation) * sin(azimuth) |
	// |V3|				  |		  -sin(elevation)         |
	Matrix& cart_from_pol(const double& magnitude, const double& azimuth, const double& elevation);

	// Returns polar from cartesian coordinates
	// magnitude = POLAR(0,0) = |V|
	// azimuth	 = POLAR(1,0) = atan2(V2, V1)
	// elevation = POLAR(2,0) = atan2(-V3, sqrt(V1^2+V2^2))
	Matrix pol_from_cart();

	// Bi-variate ellipse
	// Calculating major and minor semi-axis of ellipse and rotation angle
	//		from the symmetrical pos semi-definite MAT(2x2) matrix
	// major_semi_axis = ELLIPSE.get_loc(0,0);
	// minor_semi_axis = ELLIPSE.get_loc(1,0);
	// angle = ELLIPSE.get_loc(2,0);
	Matrix ellipse();

	// Build a square identity matrix of object 'Matrix MAT'
	Matrix& identity();

	// Returns the inverse of a square matrix MAT
	Matrix inverse();

	// Returns a 3x3 matrix row-wise from 9x1 vector
	Matrix mat33_vec9();

	// Forms matrix MAT with all elementes '1.'
	Matrix& ones();

	// Returns the skew-symmetric matrix MAT from a 3-dim vector VEC
	//		| 0 -c  b|	   |a|
	//		| c  0 -a| <-- |b|
	//		|-b  a  0|	   |c|
	Matrix skew_sym();

	// Returns the sub matrix after 'row' and 'col' have been ommitted
	Matrix sub_matrix(const int& row, const int& col);

	// Returns the transpose of a matrix Aji < Aij
	// (same as 'Matrix operator~')
	Matrix trans();

	// Returns unit vector from 3x1 vector
	Matrix univec3();

	// Returns 9x1 vector from 3x3 matrix
	Matrix vec9_mat33();

	// Forms a zero matrix MAT from object MAT(num_row, num_col)
	Matrix& zero();

	// Inequality relational operator, returns true of false
	// returns true if elements differe by more than EPS
	bool operator!=(const Matrix& B);

	// Equality relation operator
	bool operator==(const Matrix& B);

	// Scalar multiplication operator (scalar element by element multiplicaton)
	// Note: scalar must be the second operand
	Matrix operator*(const double& b);

	// Multiplication operator, returns matrix product
	Matrix operator*(const Matrix& B);

	// Scalar multiplication assignment operator
	Matrix& operator*=(const double& b);

	// Multiplication assignment operator
	Matrix& operator*=(const Matrix& B);

	// Scalar addition operator
	Matrix operator+(const double& b);

	// Addition operator, returns matrix addition
	Matrix operator+(const Matrix& B);

	// Scalar addition assignment operator
	Matrix& operator+=(const double& b);

	// Matrix addition assignment operator
	Matrix& operator+=(const Matrix& B);

	// Scalar substruction operator
	Matrix operator-(const double& b);
	
	// Substruction operator
	Matrix operator-(const Matrix& B);

	// Scalar substruction assignment operator
	Matrix& operator-=(const double& b);

	// Matrix substruction assignment operator
	Matrix& operator-=(const Matrix& B);

	// Assignement operator (deep copy)
	Matrix& operator=(const Matrix& B);

	// Returns the component(i) from vector VEC[i] or assigns a value to component VEC[i]
	double& operator[](const int& r);

	// Scalar product operator (any combination of row or column vector)
	double operator^(const Matrix& B);

	// Alternate transpose Aji < Aij
	// (same as 'Matrix trans()')
	Matrix operator~();

private:
	// Number of rows
	int m_NumRow = 0;
	// Number of columns
	int m_NumCol = 0;
	// Total number of elements
	int m_NumElem = 0;
	// Pointer to array
	double* m_Pbody = nullptr;
};

//////////////////////////////////////////////////////////
//////////////// Module utility functions ////////////////
//////////////////////////////////////////////////////////

// Returns the T.M.of the psivg -> thtvg sequence
Matrix mat2tr(const double& psivg, const double& thtvg);

// Returns the Euler T.M. of the psi -> tht -> phi sequence
Matrix mat3tr(const double& psi, const double& tht, const double& phi);

// Returns the T.M. of geographic wrt earth coordinates
Matrix cadtge(double lon, double lat);

// Returns the T.M. of earth wrt inertial coordinates
Matrix cadtei(double simulation_time);

// Returns lon, lat, alt from inertial displacement vector
Matrix cadsph(Matrix SBIE);

// Returns the T.M. of body wrt velocity coordinates (3 DoF)
// Suitable for 3 DoF back-to-turn simulations only
Matrix cadtbv(double phi, double alpha);

// Returns inertial coordinates from longitude, latitude and altitude
Matrix cadine(double lon, double lat, double alt, double time);

// Returns the sign of the function
int sign(double variable);

// Returns the angle between two 3x1 vectors
double angle(Matrix VEC1, Matrix VEC2);


// <summary>
// Table llok-up and interpolation fuynction declarations
// //////////////////////////////////////////////////////
// 
// class 'Table'
//	Stores Table data
// 
// Created on 09/07/2024 by Fedi Nabli
// </summary>
class Table
{
public:
	double* var1_values = nullptr; // Values of variable 1
	double* var2_values = nullptr; // Values of variable 2
	double* var3_values = nullptr; // Values of variable 3
	double* data = nullptr; // Table data values packaged in one array

public:
	Table() {};

	virtual ~Table()
	{
		delete var1_values;
		delete var2_values;
		delete var3_values;
		delete data;
	}

	// Getting dimension of table
	int get_dim() { return m_Dim; }

	// Getting name of table
	std::string get_name() { return m_Name; }

	// Getting 1. Independant variable dimension
	int get_var1_dim() { return m_Var1_dim; }

	// Getting 2. Independant variable dimension
	int get_var2_dim() { return m_Var2_dim; }

	// Getting 3. Independant variable dimension
	int get_var3_dim() { return m_Var3_dim; }

	// Setting dimension of table
	void set_dim(int table_dim) { m_Dim = table_dim; }

	// Setting the name of the table
	void set_name(std::string table_name) { m_Name = table_name; }

	// Setting 1. Independant variable dimension
	void set_var1_dim(int size) { m_Var1_dim = size; }

	// Setting 2. Independant variable dimension
	void set_var2_dim(int size) { m_Var2_dim = size; }

	// Setting 3. Independant variable dimension
	void set_var3_dim(int size) { m_Var3_dim = size; }

	// Setting 1. Independant variable values
	void set_var1_value(int offset, double value)
	{
		var1_values[offset] = value;
	}

	// Setting 2. Independant variable values
	void set_var2_value(int offset, double value)
	{
		var2_values[offset] = value;
	}

	// Setting 3. Independant variable values
	void set_var3_value(int offset, double value)
	{
		var3_values[offset] = value;
	}

	// Setting tabular data values
	void set_data(int offset, double value)
	{
		data[offset] = value;
	}

private:
	int m_Dim = 0; // Table dimension (1, 2 or 3)
	int m_Var1_dim = 0; // Variable 1 dimension
	int m_Var2_dim = 0; // Variable 2 dimension
	int m_Var3_dim = 0; // Variable 3 dimension

	std::string m_Name = ""; // Table name
};