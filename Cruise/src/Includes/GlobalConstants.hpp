/////////////////////////////////////////////////////////////
// FILE: 'GlobalConstants.hpp'
// Defines all global constant parameters
// 
// Created on 07/07/2024 by Fedi Nabli
/////////////////////////////////////////////////////////////

#pragma once

/* Global constants to be used in the simulation */
const double REARTH = 6370987.308;			// Mean Earth radius - m
const double WEII3 = 7.292115e-5;			// Angular rotation of Earth - Rad/s
const double RAD = 0.0174532925199432;		// Conversion factor Deg -> Rad
const double DEG = 57.2957795130823;		// Conversion factor Rad -> Deg
const double AGRAV = 9.80675445;			// Standard value of Gravity acceleration - m/s^2
const double G = 6.673e-11;					// Universal Gravitational constant - Nm^2/Kg^2
const double EARTH_MASS = 5.973e24;			// Mass of the earth - Kg
const double R = 287.053;					// Ideal gas constant (dry air) - m^2/(K*s^2)
const double PI = 3.1415927;				// Circumference of unit diameter circle
const double EPS = 1e-10;					// Machine precision error
const double SMALL = 1.e-7;					// Small Real number
const double BIG = 1e10;					// Big number
const int ILARGE = 9999;					// Large Integer number

/* Sizing of arrays */
const int CHARN = 40;						// Character numbers in variable names
const int CHARL = 150;						// Character numbers in a line
/* Verify the array sizes. If too small, dynamic memory allocations may fail! */
const int NROUND3 = 40;						// Size of 'round3' module-variable array
const int NCRUISE = 160;					// Size of 'cruise' module-variable array
const int NTARGET = 20;						// Size of 'terget' module-variable array
const int NSATELLITE = 20;					// Size of 'satellite' module-variable array
const int NEVENT = 25;						// Max number of events
const int NVAR = 15;						// Max number of variables to input at every event
