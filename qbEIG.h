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

#ifndef QBEIG_H
#define QBEIG_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "qbMatrix.h"
#include "qbVector.h"
#include "qbQR.h"

// Define error codes.
constexpr int QBEIG_MATRIXNOTSQUARE = -1;
constexpr int QBEIG_MAXITERATIONSEXCEEDED = -2;
constexpr int QBEIG_MATRIXNOTSYMMETRIC = -3;

// Function to estimate (real) eigenvalues using QR decomposition.
/* Note that this is only valid for matrices that have ALL real
	eigenvalues. The only matrices that are guaranteed to have only
	real eigenvalues are symmetric matrices. Therefore, this function
	is only guaranteed to work with symmetric matrices. */
template <typename T>
int qbEigQR(const qbMatrix2<T> &inputMatrix, std::vector<T> &eigenValues)
{
	// Make a copy of the input matrix.
	qbMatrix2<T> A = inputMatrix;

	// Verify that the input matrix is square.
	if (!A.IsSquare())
		return QBEIG_MATRIXNOTSQUARE;
		
	// Verify that the matrix is symmetric.
	if (!A.IsSymmetric())
		return QBEIG_MATRIXNOTSYMMETRIC;		
		
	// The number of eigenvalues is equal to the number of rows.
	int numRows = A.GetNumRows();
	
	// Create an identity matrix of the same dimensions.
	qbMatrix2<T> identityMatrix(numRows, numRows);
	identityMatrix.SetToIdentity();
	
	// Create matrices to store Q and R.
	qbMatrix2<T> Q (numRows, numRows);
	qbMatrix2<T> R (numRows, numRows);
	
	// Loop through each iteration.
	int maxIterations = 10e3;
	int iterationCount = 0;
	bool continueFlag = true;
	while ((iterationCount < maxIterations) && continueFlag)
	{
		// Compute the QR decomposition of A.
		int returnValue = qbQR<T>(A, Q, R);
		
		// Compute the next value of A as the product of R and Q.
		A = R * Q;
		
		/* Check if A is now close enough to being upper-triangular.
			We can do this using the IsRowEchelon() function from the 
			qbMatrix2 class. */
		if (A.IsRowEchelon())
			continueFlag = false;
						
		// Increment iterationCount.
		iterationCount++;
	}
	
	// At this point, the eigenvalues should be the diagonal elements of A.
	for (int i=0; i<numRows; ++i)
		eigenValues.push_back(A.GetElement(i,i));
	
	// Set the return status accordingly.
	if (iterationCount == maxIterations)
		return QBEIG_MAXITERATIONSEXCEEDED;
	else
		return 0;	
	
}

// Function to perform inverse power iteration method.
template <typename T>
int qbInvPIt(const qbMatrix2<T> &inputMatrix, const T &eigenValue, qbVector<T> &eigenVector)
{
	// Make a copy of the input matrix.
	qbMatrix2<T> A = inputMatrix;

	// Verify that the input matrix is square.
	if (!A.IsSquare())
		return QBEIG_MATRIXNOTSQUARE;
		
  // Setup a random number generator.
	std::random_device myRandomDevice;
  std::mt19937 myRandomGenerator(myRandomDevice());
	std::uniform_int_distribution<int> myDistribution(1.0, 10.0);
	
	/* The number of eigenvectors and eigenvalues that we will compute will be
		equal to the number of rows in the input matrix. */
	int numRows = A.GetNumRows();
	
	// Create an identity matrix of the same dimensions.
	qbMatrix2<T> identityMatrix(numRows, numRows);
	identityMatrix.SetToIdentity();
	
	// Create an initial vector, v.
	qbVector<T> v(numRows);
	for (int i=0; i<numRows; ++i)
		v.SetElement(i, static_cast<T>(myDistribution(myRandomGenerator)));
		
	// Iterate.
	int maxIterations = 100;
	int iterationCount = 0;
	T deltaThreshold = static_cast<T>(1e-9);
	T delta = static_cast<T>(1e6);
	qbVector<T> prevVector(numRows);
	qbMatrix2<T> tempMatrix(numRows, numRows);
	
	while ((iterationCount < maxIterations) && (delta > deltaThreshold))
	{
		// Store a copy of the current working vector to use for computing delta.
		prevVector = v;
		
		// Compute the next value of v.
		tempMatrix = A - (eigenValue * identityMatrix);
		tempMatrix.Inverse();
		v = tempMatrix * v;
		v.Normalize();
		
		// Compute delta.
		delta = (v - prevVector).norm();
		
		// Increment iteration count.
		iterationCount++;
	}
	
	// Return the estimated eigenvector.
	eigenVector = v;
	
	// Set the return status accordingly.
	if (iterationCount == maxIterations)
		return QBEIG_MAXITERATIONSEXCEEDED;
	else
		return 0;
		
}

// The qbEIG function (power iteration method).
template <typename T>
int qbEIG_PIt(const qbMatrix2<T> &X, T &eigenValue, qbVector<T> &eigenVector)
{
	// Make a copy of the input matrix.
	qbMatrix2<T> inputMatrix = X;

	// Verify that the input matrix is square.
	if (!inputMatrix.IsSquare())
		return QBEIG_MATRIXNOTSQUARE;
	
  // Setup a random number generator.
	std::random_device myRandomDevice;
  std::mt19937 myRandomGenerator(myRandomDevice());
	std::uniform_int_distribution<int> myDistribution(1.0, 10.0);
	
	/* The number of eigenvectors and eigenvalues that we will compute will be
		equal to the number of rows in the input matrix. */
	int numRows = inputMatrix.GetNumRows();
	
	// Create an identity matrix of the same dimensions.
	qbMatrix2<T> identityMatrix(numRows, numRows);
	identityMatrix.SetToIdentity();

	/* **************************************************************
		Compute the eigenvector.
	************************************************************** */
		
	// Create an initial vector, v.
	qbVector<T> v(numRows);
	for (int i=0; i<numRows; ++i)
		v.SetElement(i, static_cast<T>(myDistribution(myRandomGenerator)));
		
	// Loop over the required number of iterations.
	qbVector<T> v1(numRows);
	int numIterations = 1000;
	for (int i=0; i<numIterations; ++i)
	{
		v1 = inputMatrix * v;
		v1.Normalize();
		v = v1;
	}

	// Store this eigenvector.
	eigenVector = v1;
	
	/* **************************************************************
		Compute the eigenvalue corresponding to this eigenvector.
	************************************************************** */	
	
	// Compute the cumulative sum.
	T cumulativeSum = static_cast<T>(0.0);
	for (int i=1; i<numRows; ++i)
		cumulativeSum += inputMatrix.GetElement(0,i) * v1.GetElement(i);
		
	eigenValue = (cumulativeSum / v1.GetElement(0)) + inputMatrix.GetElement(0,0);

	return 0;
}

#endif
