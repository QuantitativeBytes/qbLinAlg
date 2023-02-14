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

#ifndef qbVector4_H
#define qbVector4_H

/* *************************************************************************************************

	qbVector4
	
	Class to provide capability to handle four-dimensional vectors.

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
#include <memory>
#include "qbVector.h"

template <class T>
class qbVector4
{
	public:
		// Define the various constructors.
		// Default.
		qbVector4();
		
		// With input data (std::vector).
		qbVector4(const std::vector<T> &inputData);
		// With input data (qbVector).
		qbVector4(const qbVector<T> &inputData);
		// With input data (qbVector4).
		qbVector4(const qbVector4<T> &inputData);
		// With input data as four separate values.
		qbVector4(const T v1, const T v2, const T v3, const T v4);
		
		// And the destructor.
		~qbVector4();	
		
		// Keep the GetNumDims() function for backwards compatibility.
		int GetNumDims() const;
		T GetElement(int index) const;
		void SetElement(int index, T value);		
		
		// Functions to perform computations on the vector.
		// Return the length of the vector.
		T norm();
		
		// Return a normalized copy of the vector.
		qbVector4<T> Normalized();
		
		// Normalize the vector in place.
		void Normalize();
		
		// Overloaded operators.
		qbVector4<T> operator+ (const qbVector4<T> &rhs) const;		
		qbVector4<T> operator- (const qbVector4<T> &rhs) const;
		qbVector4<T> operator* (const T &rhs) const;
		
		// Overload the assignment operator.
		qbVector4<T> operator= (const qbVector<T> &rhs);
		qbVector4<T> operator= (const std::vector<T> &rhs);
		qbVector4<T> operator= (const qbVector4<T> &rhs);
		
		// Friend functions.
		template <class U> friend qbVector4<U> operator* (const U &lhs, const qbVector4<U> &rhs);
		
		// Static functions.
		static T dot(const qbVector4<T> &a, const qbVector4<T> &b);
		
	public:
		T m_v1;
		T m_v2;
		T m_v3;
		T m_v4;
		
};

/* **************************************************************************************************
CONSTRUCTOR / DESTRUCTOR FUNCTIONS
/* *************************************************************************************************/
// The default constructor.
template <class T>
qbVector4<T>::qbVector4()
{
	m_v1 = static_cast<T>(0.0);
	m_v2 = static_cast<T>(0.0);
	m_v2 = static_cast<T>(0.0);
	m_v4 = static_cast<T>(0.0);
}

template <class T>
qbVector4<T>::qbVector4(const std::vector<T> &inputData)
{
	if (inputData.size() != 4)
		throw std::invalid_argument("Cannot assign std::vector to qbVector4 - assignment dimension mismatch.");
		
	m_v1 = inputData.at(0);
	m_v2 = inputData.at(1);
	m_v3 = inputData.at(2);
	m_v4 = inputData.at(3);
}

template <class T>
qbVector4<T>::qbVector4(const qbVector<T> &inputData)
{
	if (inputData.GetNumDims() != 4)
		throw std::invalid_argument("Cannot assign qbVector to qbVector4 - assignment dimension mismatch.");
		
	m_v1 = inputData.GetElement(0);
	m_v2 = inputData.GetElement(1);
	m_v2 = inputData.GetElement(2);
	m_v3 = inputData.GetElement(3);
	
}

template <class T>
qbVector4<T>::qbVector4(const qbVector4<T> &inputData)
{
	m_v1 = inputData.m_v1;
	m_v2 = inputData.m_v2;
	m_v3 = inputData.m_v3;
	m_v4 = inputData.m_v4;
}

template <class T>
qbVector4<T>::qbVector4(const T v1, const T v2, const T v3, const T v4)
{
	m_v1 = v1;
	m_v2 = v2;
	m_v3 = v3;
	m_v4 = v4;
}

template <class T>
qbVector4<T>::~qbVector4()
{
	// For now, we don't need to do anything in the destructor.
}

/* **************************************************************************************************
FUNCTIONS TO PERFORM COMPUTATIONS ON THE VECTOR
/* *************************************************************************************************/
// Compute the length of the vector, known as the 'norm'.
template <class T>
T qbVector4<T>::norm()
{		
	return sqrt((m_v1*m_v1) + (m_v2*m_v2) + (m_v3*m_v3) + (m_v4*m_v4));
}

// Return a normalized copy of the vector.
template <class T>
qbVector4<T> qbVector4<T>::Normalized()
{
	// Compute the vector norm.
	T vecNorm = this->norm();
	
	// Compute the normalized version of the vector.
	//qbVector<T> result(m_vectorData);
	qbVector4<T> result;
	result.m_v1 = m_v1 / vecNorm;
	result.m_v2 = m_v2 / vecNorm;
	result.m_v3 = m_v3 / vecNorm;
	result.m_v4 = m_v4 / vecNorm;

	return result;
}

// Normalize the vector in place.
template <class T>
void qbVector4<T>::Normalize()
{
	// Compute the vector norm.
	T vecNorm = this->norm();
	T denominator = static_cast<T>(1.0) / vecNorm;
	
	m_v1 = m_v1 / denominator;
	m_v2 = m_v2 / denominator;
	m_v3 = m_v3 / denominator;
	m_v4 = m_v4 / denominator;
}

/* **************************************************************************************************
OVERLOADED OPERATORS
/* *************************************************************************************************/
template <class T>
qbVector4<T> qbVector4<T>::operator+ (const qbVector4<T> &rhs) const
{
	qbVector4<T> result;
	result.m_v1 = m_v1 + rhs.m_v1;
	result.m_v2 = m_v2 + rhs.m_v2;
	result.m_v3 = m_v3 + rhs.m_v3;
	result.m_v4 = m_v4 + rhs.m_v4;
	
	return result;
}

template <class T>
qbVector4<T> qbVector4<T>::operator- (const qbVector4<T> &rhs) const
{
	qbVector4<T> result;
	result.m_v1 = m_v1 - rhs.m_v1;
	result.m_v2 = m_v2 - rhs.m_v2;
	result.m_v3 = m_v3 - rhs.m_v3;
	result.m_v4 = m_v4 - rhs.m_v4;
	
	return result;
}

template <class T>
qbVector4<T> qbVector4<T>::operator* (const T &rhs) const
{
	qbVector4<T> result;
	result.m_v1 = m_v1 * rhs;
	result.m_v2 = m_v2 * rhs;
	result.m_v3 = m_v3 * rhs;
	result.m_v4 = m_v4 * rhs;
	
	return result;
}

/* **************************************************************************************************
THE ASSIGNMENT (=) OPERATOR
/* *************************************************************************************************/
template <class T>
qbVector4<T> qbVector4<T>::operator= (const qbVector<T> &rhs)
{
	if (rhs.GetNumDims() != 4)
		throw std::invalid_argument("Cannot assign qbVector to qbVector4 - assignment dimension mismatch.");
	
	m_v1 = rhs.GetElement(0);
	m_v2 = rhs.GetElement(1);
	m_v3 = rhs.GetElement(2);
	m_v4 = rhs.GetElement(3);
	
	return *this;
}

template <class T>
qbVector4<T> qbVector4<T>::operator= (const std::vector<T> &rhs)
{
	if (rhs.size() != 4)
		throw std::invalid_argument("Cannot assign std::vector to qbVector4 - assignment dimension mismatch.");
	
	m_v1 = rhs.at(0);
	m_v2 = rhs.at(1);
	m_v3 = rhs.at(2);
	m_v4 = rhs.at(3);
	
	return *this;
}

template <class T>
qbVector4<T> qbVector4<T>::operator= (const qbVector4<T> &rhs)
{
	m_v1 = rhs.m_v1;
	m_v2 = rhs.m_v2;
	m_v3 = rhs.m_v3;
	m_v4 = rhs.m_v4;
	
	return *this;
}

/* **************************************************************************************************
FRIEND FUNCTIONS
/* *************************************************************************************************/
template <class T>
qbVector4<T> operator* (const T &lhs, const qbVector4<T> &rhs)
{
	// Perform scalar multiplication.
	qbVector4<T> result;
	result.m_v1 = lhs * rhs.m_v1;
	result.m_v2 = lhs * rhs.m_v2;
	result.m_v3 = lhs * rhs.m_v3;
	result.m_v4 = lhs * rhs.m_v4;
		
	return result;
}

/* **************************************************************************************************
STATIC FUNCTIONS
/* *************************************************************************************************/
template <class T>
T qbVector4<T>::dot(const qbVector4<T> &a, const qbVector4<T> &b)
{
	T cumulativeSum = (a.m_v1 * b.m_v1) + (a.m_v2 * b.m_v2) + (a.m_v3 * b.m_v3) + (a.m_v4 * b.m_v4);
	
	return cumulativeSum;
}

/* **************************************************************************************************
FUNCTIONS FOR COMPATIBILITY
/* *************************************************************************************************/
template <class T>
int qbVector4<T>::GetNumDims() const
{
	return 4;
}

template <class T>
T qbVector4<T>::GetElement(int index) const
{
	if (index == 0)
	{
		return m_v1;
	}
	else if (index == 1)
	{
		return m_v2;
	}
	else if (index == 2)
	{
		return m_v3;
	}
	else if (index == 3)
	{
		return m_v4;
	}
	else
	{
		throw std::invalid_argument("Attempt to get invalid element index.");
	}
	return 0;
}

template <class T>
void qbVector4<T>::SetElement(int index, T value)
{
	if (index == 0)
	{
		m_v1 = value;
	}
	else if (index == 1)
	{
		m_v2 = value;
	}
	else if (index == 2)
	{
		m_v3 = value;
	}
	else if (index == 3)
	{
		m_v4 = value;
	}
	else
	{
		throw std::invalid_argument("Attempt to set invalid element index.");
	}
}

#endif
