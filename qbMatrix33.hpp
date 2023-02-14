// This file is part of the qbLinAlg linear algebra library.

/*
MIT License
Copyright (c) 2023 Michael Bennett	

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

#ifndef qbMatrix33_H
#define qbMatrix33_H

/* *************************************************************************************************

	qbMatrix33
	
	Class to provide capability to handle two-dimensional matrices.

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
#include <exception>
#include "qbVector.h"
#include "qbVector3.hpp"
#include "qbVector4.hpp"

template <class T>
class qbMatrix33
{
	public:
		// Define the various constructors.
    qbMatrix33();
    qbMatrix33(const qbMatrix33<T> &inputMatrix);
    qbMatrix33(const std::vector<T> &inputData);

    // And the destructor.
    ~qbMatrix33();

    // Configuration methods.
    void SetToIdentity();

    // Element access methods.
    T GetElement(int row, int col) const;
    bool SetElement(int row, int col, T elementValue);
    int GetNumRows() const;
    int GetNumCols() const;
    
    // Manipulation methods.
    // Compute matrix inverse.
    bool Inverse();

    // Return the transpose.
    qbMatrix33<T> Transpose() const;
    
    // Compute determinant.
    T Determinant();

		// Overload == operator.
		bool operator== (const qbMatrix33<T>& rhs);
		bool Compare (const qbMatrix33<T>& matrix1, double tolerance);
		
		// Overload the assignment operator.
		qbMatrix33<T> operator= (const qbMatrix33<T> &rhs);

    // Overload +, - and * operators (friends).
    template <class U> friend qbMatrix33<U> operator+ (const qbMatrix33<U>& lhs, const qbMatrix33<U>& rhs);
    template <class U> friend qbMatrix33<U> operator+ (const U& lhs, const qbMatrix33<U>& rhs);
    template <class U> friend qbMatrix33<U> operator+ (const qbMatrix33<U>& lhs, const U& rhs);
    
    template <class U> friend qbMatrix33<U> operator- (const qbMatrix33<U>& lhs, const qbMatrix33<U>& rhs);
    template <class U> friend qbMatrix33<U> operator- (const U& lhs, const qbMatrix33<U>& rhs);
    template <class U> friend qbMatrix33<U> operator- (const qbMatrix33<U>& lhs, const U& rhs);
    
    template <class U> friend qbMatrix33<U> operator* (const qbMatrix33<U>& lhs, const qbMatrix33<U>& rhs);
    template <class U> friend qbMatrix33<U> operator* (const U& lhs, const qbMatrix33<U>& rhs);
    template <class U> friend qbMatrix33<U> operator* (const qbMatrix33<U>& lhs, const U& rhs);
    
    // qbMatrix33 * qbVector3.
    template <class U> friend qbVector3<U> operator* (const qbMatrix33<U>& lhs, const qbVector3<U>& rhs);      

	private:
		int Sub2Ind(int row, int col) const;
		T cofactorDeterminant(T e1, T e2, T e3, T e4);
		bool CloseEnough(T f1, T f2);

	public:
		//T *m_matrixData;
		T m_matrixData[9];
    int m_nRows = 3;
    int m_nCols = 3;
    int m_nElements = 9;
};

/* **************************************************************************************************
CONSTRUCTOR / DESTRUCTOR FUNCTIONS
/* *************************************************************************************************/
// The default constructor.
template <class T>
qbMatrix33<T>::qbMatrix33()
{

}

// The copy constructor.
template <class T>
qbMatrix33<T>::qbMatrix33(const qbMatrix33<T> &inputMatrix)
{
	for (int i=0; i<9; i++)
		m_matrixData[i] = inputMatrix.m_matrixData[i];
}

// Construct from std::vector.
template <class T>
qbMatrix33<T>::qbMatrix33(const std::vector<T> &inputData)
{
	if (inputData.size() != 9)
		throw std::invalid_argument("Cannot assign std::vector to qbMatrix33 - assignment dimension mismatch.");

	for (int i=0; i<9; ++i)
		m_matrixData[i] = inputData.at(i);
}

template <class T>
qbMatrix33<T>::~qbMatrix33()
{

}

/* **************************************************************************************************
CONFIGURATION FUNCTIONS
/* *************************************************************************************************/
// Function to convert the existing matrix into an identity matrix.
template <class T>
void qbMatrix33<T>::SetToIdentity()
{		
	for (int row=0; row<m_nRows; ++row)
	{
		for (int col=0; col<m_nCols; ++col)
		{
			if (col == row)
				m_matrixData[Sub2Ind(row,col)] = 1.0;
			else
				m_matrixData[Sub2Ind(row,col)] = 0.0;
		}
	}
}

/* **************************************************************************************************
ELEMENT FUNCTIONS
/* *************************************************************************************************/
template <class T>
T qbMatrix33<T>::GetElement(int row, int col) const
{
	int linearIndex = Sub2Ind(row, col);
	if (linearIndex >= 0)
		return m_matrixData[linearIndex];
	else
		return 0.0;

}

template <class T>
bool qbMatrix33<T>::SetElement(int row, int col, T elementValue)
{
	int linearIndex = Sub2Ind(row, col);
	if (linearIndex >= 0)
	{
		m_matrixData[linearIndex] = elementValue;
		return true;
	} 
	else 
	{
		return false;
	}
}

template <class T>
int qbMatrix33<T>::GetNumRows() const
{
	return m_nRows;
}

template <class T>
int qbMatrix33<T>::GetNumCols() const
{
	return m_nCols;
}

template <class T>
bool qbMatrix33<T>::Compare(const qbMatrix33<T>& matrix1, double tolerance)
{
	// First, check that the matrices have the same dimensions.
	int numRows1 = matrix1.m_nRows;
	int numCols1 = matrix1.m_nCols;
	if ((numRows1 != m_nRows) || (numCols1 != m_nCols))
		return false;
		
	// Loop over all the elements and compute the sum-of-squared-differences.
	double cumulativeSum = 0.0;
	for (int i=0; i<m_nElements; ++i)
	{
		T element1 = matrix1.m_matrixData[i];
		T element2 = m_matrixData[i];
		cumulativeSum += ((element1 - element2) * (element1 - element2));
	}
	double finalValue = sqrt(cumulativeSum / ((numRows1 * numCols1)-1));
	if (finalValue < tolerance)
		return true;
	else
		return false;
}

/* **************************************************************************************************
OVERLOADED OPERATOR FUNCTIONS
/* *************************************************************************************************/

/* **************************************************************************************************
THE + OPERATOR
/* *************************************************************************************************/
// matrix + matrx.
template <class T>
qbMatrix33<T> operator+ (const qbMatrix33<T>& lhs, const qbMatrix33<T>& rhs)
{
	qbMatrix33<T> result;
	for (int i=0; i<9; ++i)
		result.m_matrixData[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];

	return result;
}

// scaler + matrix
template <class T>
qbMatrix33<T> operator+ (const T& lhs, const qbMatrix33<T>& rhs)
{
	qbMatrix33<T> result;
	for (int i=0; i<9; ++i)
		result.m_matrixData[i] = lhs + rhs.m_matrixData[i];

	return result;
}

// matrix + scaler
template <class T>
qbMatrix33<T> operator+ (const qbMatrix33<T>& lhs, const T& rhs)
{
	qbMatrix33<T> result;
	for (int i=0; i<9; ++i)
		result.m_matrixData[i] = lhs.m_matrixData[i] + rhs;
			
	return result;
}

/* **************************************************************************************************
THE - OPERATOR
/* *************************************************************************************************/
// matrix - matrix
template <class T>
qbMatrix33<T> operator- (const qbMatrix33<T>& lhs, const qbMatrix33<T>& rhs)
{
	qbMatrix33<T> result;
	for (int i=0; i<9; ++i)
		result.m_matrixData[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];

	return result;   
}

// scaler - matrix
template <class T>
qbMatrix33<T> operator- (const T& lhs, const qbMatrix33<T>& rhs)
{
	qbMatrix33<T> result;
	for (int i=0; i<9; ++i)
		result.m_matrixData[i] = lhs - rhs.m_matrixData[i];

	return result;
}

// matrix - scaler
template <class T>
qbMatrix33<T> operator- (const qbMatrix33<T>& lhs, const T& rhs)
{
	qbMatrix33<T> result;
	for (int i=0; i<9; ++i)
		result.m_matrixData[i] = lhs.m_matrixData[i] - rhs;
			
	return result;
}

/* **************************************************************************************************
THE * OPERATOR
/* *************************************************************************************************/
// matrix * qbVector3
template <class T>
qbVector3<T> operator* (const qbMatrix33<T>& lhs, const qbVector3<T>& rhs)
{
	// Setup the vector for the output.
	qbVector3<T> result;
	
	// Loop over rows and columns and perform the multiplication operation element-by-element.
	for (int row=0; row<lhs.m_nRows; ++row)
	{
		T cumulativeSum = static_cast<T>(0.0);
		cumulativeSum += (lhs.GetElement(row,0) * rhs.m_x);
		cumulativeSum += (lhs.GetElement(row,1) * rhs.m_y);
		cumulativeSum += (lhs.GetElement(row,2) * rhs.m_z);

		result.SetElement(row, cumulativeSum);
	}
	
	return result;
}

// scaler * matrix
template <class T>
qbMatrix33<T> operator* (const T& lhs, const qbMatrix33<T>& rhs)
{
	qbMatrix33<T> result;
	for (int i=0; i<9; ++i)
		result.m_matrixData[i] = lhs * rhs.m_matrixData[i];
			
	return result;
}

// matrix * scaler
template <class T>
qbMatrix33<T> operator* (const qbMatrix33<T>& lhs, const T& rhs)
{
	qbMatrix33<T> result;
	for (int i=0; i<9; ++i)
		result.m_matrixData[i] = lhs.m_matrixData[i] * rhs;
			
	return result;
}

// matrix * matrix
template <class T>
qbMatrix33<T> operator* (const qbMatrix33<T>& lhs, const qbMatrix33<T>& rhs)
{
	int r_numRows = rhs.m_nRows;
	int r_numCols = rhs.m_nCols;
	int l_numRows = lhs.m_nRows;
	int l_numCols = lhs.m_nCols;

	// This is the standard matrix multiplication condition.
	// The output will be the same size as the RHS.
	//T *tempResult = new T[lhs.m_nRows * rhs.m_nCols];
	qbMatrix33<T> result;

	// Loop through each row of the LHS.
	for (int lhsRow=0; lhsRow<l_numRows; lhsRow++)
	{
		// Loop through each column on the RHS.
		for (int rhsCol=0; rhsCol<r_numCols; rhsCol++)
		{
			T elementResult = static_cast<T>(0.0);
			
			// Loop through each element of this LHS row.
			for (int lhsCol=0; lhsCol<l_numCols; lhsCol++)
			{
				// Compute the LHS linear index.
				int lhsLinearIndex = (lhsRow * l_numCols) + lhsCol;
				
				// Compute the RHS linear index (based on LHS col).
				// rhs row number equal to lhs column number.
				int rhsLinearIndex = (lhsCol * r_numCols) + rhsCol;
				
				// Perform the calculation on these elements.
				elementResult += (lhs.m_matrixData[lhsLinearIndex] * rhs.m_matrixData[rhsLinearIndex]);
			}
			// Store the result.
			int resultLinearIndex = (lhsRow * r_numCols) + rhsCol;
			result.m_matrixData[resultLinearIndex] = elementResult;
			//tempResult[resultLinearIndex] = elementResult;
		}		
	}
	return result;
}

/* **************************************************************************************************
THE == OPERATOR
/* *************************************************************************************************/
template <class T>
bool qbMatrix33<T>::operator== (const qbMatrix33<T>& rhs)
{	
	// Check if the elements are equal.
	bool flag = true;
	for (int i=0; i<9; ++i)
	{
		if (!CloseEnough(this->m_matrixData[i], rhs.m_matrixData[i]))
			flag = false;
	}
	return flag;
}

/* **************************************************************************************************
THE ASSIGNMENT (=) OPERATOR
/* *************************************************************************************************/
template <class T>
qbMatrix33<T> qbMatrix33<T>::operator= (const qbMatrix33<T> &rhs)
{
	// Make sure we're not assigning to ourself.
	if (this != &rhs)
	{
		for (int i=0; i<m_nElements; ++i)
			m_matrixData[i] = rhs.m_matrixData[i];
	}
	
	return *this;
}

/* **************************************************************************************************
TRANSPOSE
/* *************************************************************************************************/
template <class T>
qbMatrix33<T> qbMatrix33<T>::Transpose() const
{
	qbMatrix33<double> result;

	for (int i=0; i<3; ++i)
	{
		for (int j=0; j<3; ++j)
		{
			result.SetElement(j, i, this->GetElement(i, j));
		}
	}	

	return result;
}

/* **************************************************************************************************
INVERSE
/* *************************************************************************************************/
template <class T>
bool qbMatrix33<T>::Inverse()
{
	// Compute the adjucate matrix.
	qbMatrix33<T> adjugate;
	adjugate.SetElement(0, 0, cofactorDeterminant(m_matrixData[4], m_matrixData[5], m_matrixData[7], m_matrixData[8]));
	adjugate.SetElement(0, 1, -cofactorDeterminant(m_matrixData[3], m_matrixData[5], m_matrixData[6], m_matrixData[8]));
	adjugate.SetElement(0, 2, cofactorDeterminant(m_matrixData[3], m_matrixData[4], m_matrixData[6], m_matrixData[7]));
	
	adjugate.SetElement(1, 0, -cofactorDeterminant(m_matrixData[1], m_matrixData[2], m_matrixData[7], m_matrixData[8]));
	adjugate.SetElement(1, 1, cofactorDeterminant(m_matrixData[0], m_matrixData[2], m_matrixData[6], m_matrixData[8]));
	adjugate.SetElement(1, 2, -cofactorDeterminant(m_matrixData[0], m_matrixData[1], m_matrixData[6], m_matrixData[7]));
	
	adjugate.SetElement(2, 0, cofactorDeterminant(m_matrixData[1], m_matrixData[2], m_matrixData[4], m_matrixData[5]));
	adjugate.SetElement(2, 1, -cofactorDeterminant(m_matrixData[0], m_matrixData[2], m_matrixData[3],m_matrixData[5]));
	adjugate.SetElement(2, 2, cofactorDeterminant(m_matrixData[0], m_matrixData[1], m_matrixData[3], m_matrixData[4]));
	
	// And transpose.
	qbMatrix33<T> adjT = adjugate.Transpose();
	
	// Compute the inverse.
	T determinant = Determinant();
	if (CloseEnough(determinant, 0.0))
		return false;
	
	qbMatrix33<T> result = (1.0 / determinant) * adjT;
	
	// And store 'in place'.
	for (int i=0; i<9; ++i)
		m_matrixData[i] = result.m_matrixData[i];
	
	return true;
}

/* **************************************************************************************************
DETERMINANT
/* *************************************************************************************************/
template <class T>
T qbMatrix33<T>::Determinant()
{
	T result =	(m_matrixData[0] * m_matrixData[4] * m_matrixData[8]) + 
							(m_matrixData[1] * m_matrixData[5] * m_matrixData[6]) +
							(m_matrixData[2] * m_matrixData[3] * m_matrixData[7]) -
							(m_matrixData[2] * m_matrixData[4] * m_matrixData[6]) -
							(m_matrixData[1] * m_matrixData[3] * m_matrixData[8]) -
							(m_matrixData[0] * m_matrixData[5] * m_matrixData[7]);
	
	return result;
}

/* **************************************************************************************************
PRIVATE FUNCTIONS
/* *************************************************************************************************/
// Function to return the linear index corresponding to the supplied row and column values.
template <class T>
int qbMatrix33<T>::Sub2Ind(int row, int col) const
{
	if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
		return (row * m_nCols) + col;
	else
		return -1;
}

// Function to compute the determinant of a 2x2 cofactor matrix.
template <class T>
T qbMatrix33<T>::cofactorDeterminant(T e1, T e2, T e3, T e4)
{
	return (e1 * e4) - (e2 * e3);
}

template <class T>
bool qbMatrix33<T>::CloseEnough(T f1, T f2)
{
    return fabs(f1-f2) < 1e-9;
}

#endif
