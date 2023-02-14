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

#ifndef QBVECTOR_H
#define QBVECTOR_H

/* *************************************************************************************************

	qbVector
	
	Class to provide capability to handle vectors.

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

template <class T>
class qbVector
{
	public:
		// Define the various constructors.
		// Default.
		qbVector();
		
		// With a single integer specifying the number of dimensions.
		qbVector(int numDims);
		
		// With input data (std::vector).
		qbVector(std::vector<T> inputData);
		
		// And the destructor.
		~qbVector();	
		
		// Functions to return parameters of the vector.
		int GetNumDims() const;
		
		// Functions to handle elements of the vector.
		T GetElement(int index) const;
		void SetElement(int index, T value);
		
		// Functions to perform computations on the vector.
		// Return the length of the vector.
		T norm();
		
		// Return a normalized copy of the vector.
		qbVector<T> Normalized();
		
		// Normalize the vector in place.
		void Normalize();
		
		// Overloaded operators.
		qbVector<T> operator+ (const qbVector<T> &rhs) const;		
		qbVector<T> operator- (const qbVector<T> &rhs) const;
		qbVector<T> operator* (const T &rhs) const;
		
		// Friend functions.
		template <class U> friend qbVector<U> operator* (const U &lhs, const qbVector<U> &rhs);
		
		// Static functions.
		static T dot(const qbVector<T> &a, const qbVector<T> &b);
		static qbVector<T> cross(const qbVector<T> &a, const qbVector<T> &b);
		
	private:
		std::vector<T> m_vectorData;
		int m_nDims;
		
};

/* **************************************************************************************************
CONSTRUCTOR / DESTRUCTOR FUNCTIONS
/* *************************************************************************************************/
// The default constructor.
template <class T>
qbVector<T>::qbVector()
{
	m_nDims = 0;
	m_vectorData = std::vector<T>();
}

template <class T>
qbVector<T>::qbVector(int numDims)
{
	m_nDims = numDims;
	m_vectorData = std::vector<T>(numDims, static_cast<T>(0.0));
}

template <class T>
qbVector<T>::qbVector(std::vector<T> inputData)
{
	m_nDims = inputData.size();
	m_vectorData = inputData;
}

template <class T>
qbVector<T>::~qbVector()
{
	// For now, we don't need to do anything in the destructor.
}

/* **************************************************************************************************
FUNCTIONS TO RETURN PARAMETERS
/* *************************************************************************************************/
template <class T>
int qbVector<T>::GetNumDims() const
{
	return m_nDims;
}

/* **************************************************************************************************
FUNCTIONS TO HANDLE ELEMENTS OF THE VECTOR
/* *************************************************************************************************/
template <class T>
T qbVector<T>::GetElement(int index) const
{
	return m_vectorData[index];
}

template <class T>
void qbVector<T>::SetElement(int index, T value)
{
	m_vectorData[index] = value;
}

/* **************************************************************************************************
FUNCTIONS TO PERFORM COMPUTATIONS ON THE VECTOR
/* *************************************************************************************************/
// Compute the length of the vector,known as the 'norm'.
template <class T>
T qbVector<T>::norm()
{
	T cumulativeSum = static_cast<T>(0.0);
	for (int i=0; i<m_nDims; ++i)
		cumulativeSum += (m_vectorData.at(i) * m_vectorData.at(i));
		
	return sqrt(cumulativeSum);
}

// Return a normalized copy of the vector.
template <class T>
qbVector<T> qbVector<T>::Normalized()
{
	// Compute the vector norm.
	T vecNorm = this->norm();
	
	// Compute the normalized version of the vector.
	qbVector<T> result(m_vectorData);
	return result * (static_cast<T>(1.0) / vecNorm);
}

// Normalize the vector in place.
template <class T>
void qbVector<T>::Normalize()
{
	// Compute the vector norm.
	T vecNorm = this->norm();
	
	// Compute the elements of the normalized version of the vector.
	for (int i=0; i<m_nDims; ++i)
	{
		T temp = m_vectorData.at(i) * (static_cast<T>(1.0) / vecNorm);
		m_vectorData.at(i) = temp;
	}
}

/* **************************************************************************************************
OVERLOADED OPERATORS
Note the changes that have been made since the videos about this first code were made. The original
code has been left, but commented out. Ultimately it was necessary to make some changes to 
improve the performance of these functions in terms of run time.
See this episode of my ray tracing in C++ series for further information:
https://youtu.be/-5kLk7_bs0U
/* *************************************************************************************************/
template <class T>
qbVector<T> qbVector<T>::operator+ (const qbVector<T> &rhs) const
{
	// Check that the number of dimensions match.
	if (m_nDims != rhs.m_nDims)
		throw std::invalid_argument("Vector dimensions do not match.");
	
	/*
	std::vector<T> resultData;
	for (int i=0; i<m_nDims; ++i)
		resultData.push_back(m_vectorData.at(i) + rhs.m_vectorData.at(i));
		
	qbVector<T> result(resultData);
	return result;
	*/
	qbVector<T> resultData (m_nDims);
	for (int i=0; i<m_nDims; ++i)
		resultData.SetElement(i, (m_vectorData.at(i) + rhs.m_vectorData.at(i)));
		
	return resultData;	
}

template <class T>
qbVector<T> qbVector<T>::operator- (const qbVector<T> &rhs) const
{
	// Check that the number of dimensions match.
	if (m_nDims != rhs.m_nDims)
		throw std::invalid_argument("Vector dimensions do not match.");
	
	/*
	std::vector<T> resultData;
	for (int i=0; i<m_nDims; ++i)
		resultData.push_back(m_vectorData.at(i) - rhs.m_vectorData.at(i));
		
	qbVector<T> result(resultData);
	return result;
	*/
	qbVector<T> resultData (m_nDims);
	for (int i=0; i<m_nDims; ++i)
		resultData.SetElement(i, (m_vectorData.at(i) - rhs.m_vectorData.at(i)));
		
	return resultData;	
}

template <class T>
qbVector<T> qbVector<T>::operator* (const T &rhs) const
{
	// Perform scalar multiplication.
	/*
	std::vector<T> resultData;
	for (int i=0; i<m_nDims; ++i)
		resultData.push_back(m_vectorData.at(i) * rhs);
		
	qbVector<T> result(resultData);
	return result;
	*/
	qbVector<T> resultData (m_nDims);
	for (int i=0; i<m_nDims; ++i)
		resultData.SetElement(i, (m_vectorData.at(i) * rhs));
		
	return resultData;	
}

/* **************************************************************************************************
FRIEND FUNCTIONS
/* *************************************************************************************************/
template <class T>
qbVector<T> operator* (const T &lhs, const qbVector<T> &rhs)
{
	// Perform scalar multiplication.
	/*
	std::vector<T> resultData;
	for (int i=0; i<rhs.m_nDims; ++i)
		resultData.push_back(lhs * rhs.m_vectorData.at(i));
		
	qbVector<T> result(resultData);
	return result;
	*/
	std::vector<T> resultData (rhs.m_nDims);
	for (int i=0; i<rhs.m_nDims; ++i)
		resultData.at(i) = (lhs * rhs.m_vectorData.at(i));
		
	qbVector<T> result(resultData);
	return result;	
}

/* **************************************************************************************************
STATIC FUNCTIONS
/* *************************************************************************************************/
template <class T>
T qbVector<T>::dot(const qbVector<T> &a, const qbVector<T> &b)
{
	// Check that the number of dimensions match.
	if (a.m_nDims != b.m_nDims)
		throw std::invalid_argument("Vector dimensions must match for the dot-product to be computed.");
	
	// Compute the dot product.
	T cumulativeSum = static_cast<T>(0.0);
	for (int i=0; i<a.m_nDims; ++i)
		cumulativeSum += a.m_vectorData.at(i) * b.m_vectorData.at(i);
	
	return cumulativeSum;
}

template <class T>
qbVector<T> qbVector<T>::cross(const qbVector<T> &a, const qbVector<T> &b)
{
	// Check that the number of dimensions match.
	if (a.m_nDims != b.m_nDims)
		throw std::invalid_argument("Vector dimensions must match for the cross-product to be computed.");
	
	// Check that the number of dimensions is 3.
	/* Although the cross-product is also defined for 7 dimensions, we are
		not going to consider that case at this time. */
	if (a.m_nDims != 3)
		throw std::invalid_argument("The cross-product can only be computed for three-dimensional vectors.");
	
	// Compute the cross product.
	std::vector<T> resultData;
	resultData.push_back((a.m_vectorData.at(1) * b.m_vectorData.at(2)) - (a.m_vectorData.at(2) * b.m_vectorData.at(1)));
	resultData.push_back(-((a.m_vectorData.at(0) * b.m_vectorData.at(2)) - (a.m_vectorData.at(2) * b.m_vectorData.at(0))));
	resultData.push_back((a.m_vectorData.at(0) * b.m_vectorData.at(1)) - (a.m_vectorData.at(1) * b.m_vectorData.at(0)));

	qbVector<T> result(resultData);
	return result;
}

#endif
