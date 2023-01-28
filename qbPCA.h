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

#ifndef QBPCA_H
#define QBPCA_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <algorithm>

#include "qbMatrix.h"
#include "qbVector.h"
#include "qbEIG.h"

// Define error codes.
constexpr int QBPCA_MATRIXNOTSQUARE = -1;
constexpr int QBPCA_MATRIXNOTSYMMETRIC = -2;

namespace qbPCA
{

// Function to compute the column means.
template <typename T>
std::vector<T> ComputeColumnMeans(const qbMatrix2<T> &inputData)
{
	// Determine the size of the input data.
	int numRows = inputData.GetNumRows();
	int numCols = inputData.GetNumCols();
	
	// Create a vector for output.
	std::vector<T> output;	
	
	// Loop through and compute means.
	for (int j=0; j<numCols; ++j)
	{
		T cumulativeSum = static_cast<T>(0.0);
		for (int i=0; i<numRows; ++i)
			cumulativeSum += inputData.GetElement(i,j);
			
		output.push_back(cumulativeSum / static_cast<T>(numRows));
	}
	
	return output;
}

// Function to subtract the column means.
template <typename T>
void SubtractColumnMeans(qbMatrix2<T> &inputData, std::vector<T> &columnMeans)
{
	// Determine the size of the input data.
	int numRows = inputData.GetNumRows();
	int numCols = inputData.GetNumCols();
	
	// Loop through and subtract the means.
	for (int j=0; j<numCols; ++j)
	{
		for (int i=0; i<numRows; ++i)
			inputData.SetElement(i,j, inputData.GetElement(i,j) - columnMeans.at(j));
	}	
}

// Function to compute the covaraince matrix.
template <typename T>
qbMatrix2<T> ComputeCovariance(const qbMatrix2<T> &X)
{
	/* Compute the covariance matrix.
		Note that here we use X'X, rather than XX' as is the usual case.
		This is because we are requiring our data to be arranged with one 
		column (p) for each variable, with one row (k) for each observation. If
		we computed XX', the result would be a [k x k] matrix. The covariance
		matrix should be [p x p], so we need to transpose, hence the use of
		X'X. */
	int numRows = X.GetNumRows();
	qbMatrix2<T> covX = (static_cast<T>(1.0) / static_cast<T>(numRows - 1)) * (X.Transpose() * X);
	return covX;
}

// Function to compute the eigenvectors of the covariance matrix.
template <typename T>
int ComputeEigenvectors(const qbMatrix2<T> &covarianceMatrix, qbMatrix2<T> &eigenvectors)
{
	// Copy the input matrix.
	qbMatrix2<T> X = covarianceMatrix;

	// The covariance matrix must be square and symmetric.
	if (!X.IsSquare())
		return QBPCA_MATRIXNOTSQUARE;
		
	// Verify that the matrix is symmetric.
	if (!X.IsSymmetric())
		return QBPCA_MATRIXNOTSYMMETRIC;
		
	// Compute the eignvalues.
	std::vector<T> eigenValues;
	int returnStatus = qbEigQR(X, eigenValues);

	// Sort the eigenvalues.
	std::sort(eigenValues.begin(), eigenValues.end());
	std::reverse(eigenValues.begin(), eigenValues.end());

	// Compute the eigenvector for each eigenvalue.
	qbVector<T> eV(X.GetNumCols());
	qbMatrix2<T> eVM(X.GetNumRows(), X.GetNumCols());
	for (int j=0; j<eigenValues.size(); ++j)
	{
		T eig = eigenValues.at(j);
		int returnStatus2 = qbInvPIt<T>(X, eig, eV);
		for (int i=0; i<eV.GetNumDims(); ++i)
			eVM.SetElement(i, j, eV.GetElement(i));
	}
	
	// Return the eigenvectors.
	eigenvectors = eVM;

	// Return the final return status.	
	return returnStatus;
}

/* Function to compute the principal components of the supplied data. */
template <typename T>
int qbPCA(const qbMatrix2<T> &inputData, qbMatrix2<T> &outputComponents)
{
	// Make a copy of the input matrix.
	qbMatrix2<T> X = inputData;
	
	// Compute the mean of each column of X.
	std::vector<T> columnMeans = ComputeColumnMeans(X);
	
	// Subtract the column means from the data.
	SubtractColumnMeans<T>(X, columnMeans);
	
	// Compute the covariance matrix.
	qbMatrix2<T> covX = ComputeCovariance(X);
	
	// Compute the eigenvectors.
	qbMatrix2<T> eigenvectors;
	int returnStatus = ComputeEigenvectors(covX, eigenvectors);
	
	// Return the output.
	outputComponents = eigenvectors;
	
	return returnStatus;
}

}

#endif
