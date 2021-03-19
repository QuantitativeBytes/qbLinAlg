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
constexpr int QBEIG_MAXITERATIONSEXCEEDED = -2;

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
	T deltaThreshold = static_cast<T>(1e-3);
	T delta = static_cast<T>(1e6);
	qbVector<T> prevVector(numRows);
	
	while ((iterationCount < maxIterations) && (delta > deltaThreshold))
	{
		// Store a copy of the current working vector to use for computing delta.
		prevVector = v;
		
		// Compute the next value of v.
		qbMatrix2<T> tempMatrix = A - (eigenValue * identityMatrix);
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
