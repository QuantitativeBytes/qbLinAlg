// This file is part of the qbLinAlg linear algebra library.
// Copyright (c) 2021 Michael Bennett
// MIT license

#ifndef QBLINSOLVE_H
#define QBLINESOLVE_H

/* *************************************************************************************************

	qbLinSolve
	
	Function to solve a system of linear equations in the form of y = X*beta, where we
	want to solve for beta.
	
	*** INPUTS ***
	
	aMatrix		qbMatrix2<T>		The matrix of independent variables (X in the above equation).
	bVector		qbVector<T>		The vector of dependent variables (y in the above equation).
	resultVec	qbVector<T>		The vector of unknown parameters (beta in the above equation).
						The final solution is returned in this vector.
															
	*** OUTPUTS ***
	
	INT				Flag indicating success or failure of the process.
						1 Indicates success.
						-1 indicates failure due to there being no unique solution (infinite solutions).
						-2 indicates failure due to there being no solution.
								
	Uses Gaussian elimination on the augmented matrix, followed by back substitution.

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

#include "qbMatrix.h"
#include "qbVector.h"

// Define error codes.
constexpr int QBLINSOLVE_NOUNIQUESOLUTION = -1;
constexpr int QBLINSOLVE_NOSOLUTIONS = -2;

// The qbLinSolve function.
template <typename T>
int qbLinSolve(const qbMatrix2<T> &aMatrix, const qbVector<T> &bVector, qbVector<T> &resultVec)
{
	// Make a copy of the input matrix, aMatrix.
	// We will use this to create the augmented matrix, so we have
	// to make a copy.
	qbMatrix2<T> inputMatrix = aMatrix;
	
	// Compute the rank of the original matrix.
	int originalRank = inputMatrix.Rank();

	/* Combine inputMatrix and bVector together into a single matrix,
		ready for using Gaussian elimination to reduce to 
		row-echelon form. */
	
	// Extract data from bVector.
	int numDims = bVector.GetNumDims();
	std::vector<T> bVecData;
	for (int i=0; i<numDims; ++i)
		bVecData.push_back(bVector.GetElement(i));
		
	// Use this to create a qbMatrix2 object with the same data (nx1).
	qbMatrix2<T> bMatrix(numDims, 1, bVecData);
	
	// Combine the two matrices together.
	inputMatrix.Join(bMatrix);
	
	/* Use Gaussian elmination to convert to row-echelon form. */
	qbMatrix2<T> rowEchelonMatrix = inputMatrix.RowEchelon();
	
	/* Comute the rank of the augmented matrix.
		Note that we do this after performing Gaussian elimination to
		reduce the matrix to row echelon form so that if this was 
		successful, there is no need to repeat this operation twice. */
	int augmentedRank = rowEchelonMatrix.Rank();
	
	/* ********************************************************************* 
		Test the two ranks to determine the nature of the system we
		are dealing with. The conditions are as follows:
		
		n = number of rows.
		
		1) originalRank = augmentedRank = n	=> A unique solution exists.
		2) originalRank = augmentedRank < n	=> An infinite number of solutions exist.
		3) originalRank < augmentedRank			=> No solutions exist.  
		********************************************************************* */
	if ((originalRank == augmentedRank) && (originalRank < inputMatrix.GetNumRows()))
	{
		return QBLINSOLVE_NOUNIQUESOLUTION;
	}
	else if (originalRank < augmentedRank)
	{
		return QBLINSOLVE_NOSOLUTIONS;
	}
	else
	{
		/* Create a qbVector object to store the output. Initially we will
			populate this with the data from bVecData, but we are going to modify
			the elements as we compute them. */
		qbVector<T> output(bVecData);
		
		// Now use back-substitution to compute the result.
		int numRows = rowEchelonMatrix.GetNumRows();
		int numCols = rowEchelonMatrix.GetNumCols();
		int startRow = numRows-1;
		
		// Loop over the rows, in reverse order.
		for (int i=startRow; i>=0; --i)
		{
			// Extract the currentResult for this row.
			T currentResult = rowEchelonMatrix.GetElement(i, numCols-1);
	
			// Compute the cumulative sum.
			T cumulativeSum = static_cast<T>(0.0);
			for (int j=i+1; j<numRows; ++j)
			{
				cumulativeSum += (rowEchelonMatrix.GetElement(i,j) * output.GetElement(j));
			}
			
			// Compute the answer.
			T finalAnswer = (currentResult - cumulativeSum) / rowEchelonMatrix.GetElement(i,i);
			
			// And store.
			output.SetElement(i, finalAnswer);
			
		}
		
		// Return the output.
		resultVec = output;
		
	}
		
	return 1;	
}

#endif
