#ifndef QBEIG_H
#define QBEIG_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "qbMatrix.h"
#include "qbVector.h"

// Define error codes.
constexpr int QBEIG_MATRIXNOTSQUARE = -1;

// The qbEIG function.
template <typename T>
int qbEIG_PIt(const qbMatrix2<T> X, T &eigenValue, qbVector<T> &eigenVector)
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
