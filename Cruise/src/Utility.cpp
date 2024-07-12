/////////////////////////////////////////////////////////////
// FILE: 'Utility.cpp'
// 
// Created on 02/07/2024 by Fedi Nabli
/////////////////////////////////////////////////////////////

#include "Includes/Utility.hpp"

////////////////////////////////////////////////////////////////////////////
////////////////// Table look-up and Integration functions /////////////////
////////////////////////////////////////////////////////////////////////////

// Single independant variable look-up
// Constant extrapolation at the upper end, slope extrapolation at the lower end
double Datadeck::look_up(std::string name, double value1)
{
	// Finding slot of table in table pointer array (Table** m_TablePtr)
	int slot(-1);
	std::string tbl_name;

	do
	{
		slot++;
		tbl_name = get_tbl(slot)->get_name();
	} while (name != tbl_name);

	// Getting table index locator of discrete value just below of variable value
	int var1_dim = get_tbl(slot)->get_var1_dim();
	int loc1 = find_index(var1_dim - 1, value1, get_tbl(slot)->var1_values);

	// Using max discrete value if value is outside table
	if (loc1 == (var1_dim - 1))
		return get_tbl(slot)->data[loc1];

	return interpolate(loc1, loc1 + 1, slot, value1);
}

// Two independant variables look-up
// Constant extrapolation at the upper end, slope extrapolation at the lower end
double Datadeck::look_up(std::string name, double value1, double value2)
{
	// Finding slot of table in table pointer array (Table** m_TablePtr)
	int slot(-1);
	std::string tbl_name;

	do
	{
		slot++;
		tbl_name = get_tbl(slot)->get_name();
	} while (name != tbl_name);

	// Getting table index (offset) locator of discrete value just below or equal of the variable value
	int var1_dim = get_tbl(slot)->get_var1_dim();
	int loc1 = find_index(var1_dim - 1, value1, get_tbl(slot)->var1_values);

	int var2_dim = get_tbl(slot)->get_var2_dim();
	int loc2 = find_index(var2_dim - 1, value2, get_tbl(slot)->var2_values);

	return interpolate(loc1, loc1 + 1, loc2, loc2 + 1, slot, value1, value2);
}

// Three independant variables look-up
// Constant extrapolation at the upper end, slope extrapolation at the lower end
double Datadeck::look_up(std::string name, double value1, double value2, double value3)
{
	// Finding slot of table in table pointer array (Table** m_TablePtr)
	int slot(-1);
	std::string tbl_name;

	do
	{
		slot++;
		tbl_name = get_tbl(slot)->get_name();
	} while (name != tbl_name);

	// Getting table index (offset) locator of discrete value just below of variable value
	int var1_dim = get_tbl(slot)->get_var1_dim();
	int loc1 = find_index(var1_dim - 1, value1, get_tbl(slot)->var1_values);

	int var2_dim = get_tbl(slot)->get_var2_dim();
	int loc2 = find_index(var2_dim - 1, value2, get_tbl(slot)->var2_values);

	int var3_dim = get_tbl(slot)->get_var3_dim();
	int loc3 = find_index(var3_dim - 1, value3, get_tbl(slot)->var3_values);

	return interpolate(loc1, loc1 + 1, loc2, loc2 + 1, loc3, loc3 + 1, slot, value1, value2, value3);
}

// Table index finder
// This is a binary search method it is O(logN)
// * Returns array locator of the discrete_variable just below variable
// * Keeps max or min array locator if variable is outside those max or min
int Datadeck::find_index(int max, double value, double* list)
{
	if (value >= list[max])
		return max;

	else if (value <= list[0])
		return 0;

	else
	{
		int index = 0;
		int mid;

		while (index <= max)
		{
			mid = (index + max) / 2;
			if (value < list[mid])
				max = mid - 1;
			else if (value > list[mid])
				index = mid + 1;
			else
				return mid;
		}

		return max;
	}
}

///////////////////////////////////////////////////////////////////////////////
// Linear one-dimensional interpolation
// Data deck must contain table in the following format:
//
//	X1       Table Values(X1)
//
//	X11		Y11
//	X12		Y12
//	X13		Y13
//           
//	Constant extrapolation beyond max values of X1
//	Slope extrapolation beyond min values of X1
///////////////////////////////////////////////////////////////////////////////
double Datadeck::interpolate(int ind1, int ind2, int slot, double val)
{
	double dx(0), dy(0);
	double dumx(0);

	double diff = val - get_tbl(slot)->var1_values[ind1];
	dx = get_tbl(slot)->var1_values[ind2] - get_tbl(slot)->var1_values[ind1];
	dy = get_tbl(slot)->data[ind2] - get_tbl(slot)->data[ind1];

	if (dx > EPS)
		dumx = diff / dx;
	dy = dumx * dy;

	return get_tbl(slot)->data[ind1] + dy;
}

///////////////////////////////////////////////////////////////////////////////
// Linear, two-dimensional interpolation
// File must contain table in the following form:
//
//  X1  X2  //Table Values(X1-row, X2-column)
//            ---------------
//  X11 X21   |Y11  Y12  Y13| 
//  X12 X22   |Y21  Y22  Y23|    <- data
//  X13 X23   |Y31  Y32  Y33| 
//  X14       |Y41  Y42  Y43| 
//            ---------------
//	Constant extrapolation beyond max values of X1 and X2
//	Slope extrapolation beyond min values of X1 and X2
///////////////////////////////////////////////////////////////////////////////
double Datadeck::interpolate(int ind10, int ind11, int ind20, int ind21, int slot, double value1, double value2)
{
	double dx1(0), dx2(0);
	double dumx1(0), dumx2(0);

	int var1_dim = get_tbl(slot)->get_var1_dim();
	int var2_dim = get_tbl(slot)->get_var2_dim();

	double diff1 = value1 - get_tbl(slot)->var1_values[ind10];
	double diff2 = value2 - get_tbl(slot)->var2_values[ind20];

	if (ind10 == (var1_dim - 1)) // Assures constant upper extrapolation of first variable
		ind11 = ind10;
	else
		dx1 = get_tbl(slot)->var1_values[ind11] - get_tbl(slot)->var1_values[ind10];

	if (ind20 == (var2_dim - 1)) // Assures constant upper extrapolation of second variable
		ind21 = ind20;
	else
		dx2 = get_tbl(slot)->var2_values[ind21] - get_tbl(slot)->var2_values[ind20];

	if (dx1 > EPS)
		dumx1 = diff1 / dx1;
	if (dx2 > EPS)
		dumx2 = diff2 / dx2;

	double y11 = get_tbl(slot)->data[ind10 * var2_dim + ind20];
	double y12 = get_tbl(slot)->data[ind10 * var2_dim + ind21];
	double y21 = get_tbl(slot)->data[ind11 * var2_dim + ind20];
	double y22 = get_tbl(slot)->data[ind11 * var2_dim + ind21];
	double y1 = dumx1 * (y21 - y11) + y11;
	double y2 = dumx2 * (y22 - y11) + y12;

	return dumx2 * (y2 - y1) + y1;
}


