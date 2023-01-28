// This file is part of the qbLinAlg linear algebra library.
/*
MIT License
Copyright (c) 2021 Michael Bennett

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, 
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or 
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef QBLSQ_H
#define QBLSQ_H

/* *************************************************************************************************

	qbLSQ
	
	Function to solve a system of linear equations using a least squares approach to handle systems
	where there are more equations (observations) than unknowns. Assumes that the system is in the
	form of y = X*beta.
	
	*** INPUTS ***
	
	Xin		qbMatrix2<T>		The matrix of independent variables (X in the above equation).
	yin		qbVector<T>		The vector of dependent variables (y in the above equation).
	result		qbVector<T>		The vector of unknown parameters (beta in the above equation).
						The final solution is returned in this vector.
															
	*** OUTPUTS ***
	
	INT				Flag indicating success or failure of the process.
						1 Indicates success.
						-1 indicates failure due to there being no computable inverse.

	Created as part of the qbLinAlg linear algebra library, which is intended to be primarily for
	educational purposes. For more details, see the corresponding videos on the QuantitativeBytes
	YouTube channel at:
	
	www.youtube.com/c/QuantitativeBytes								

	************************************************************************************************* */

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "qbVector.h"
#include "qbMatrix.h"

// Define error codes.
constexpr int QBLSQ_NOINVERSE = -1;

// The qbLSQ function.
template <typename T>
int qbLSQ(const qbMatrix2<T> &Xin, const qbVector<T> &yin, qbVector<T> &result)
{
	// Firstly, make a copy of X and y.
	qbMatrix2<T> X = Xin;
	qbVector<T> y = yin;
	
	// Compute the tranpose of X.
	qbMatrix2<T> XT = X.Transpose();
	
	// Compute XTX.
	qbMatrix2<T> XTX = XT * X;
	
	// Compute the inverse of this.
	if (!XTX.Inverse())
	{
		// We were unable to compute the inverse.
		return QBLSQ_NOINVERSE;
	}
	
	// Multiply the inverse by XT.
	qbMatrix2<T> XTXXT = XTX * XT;
	
	// And multiply by y to get the final result.
	result = XTXXT * y;
	
	return 1;
}

#endif
