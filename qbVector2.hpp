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

#ifndef QBVECTOR2_H
#define QBVECTOR2_H

/* *************************************************************************************************

	qbVector2
	
	Class to provide capability to handle two-dimensional vectors.

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
class qbVector2
{
	public:
		// Define the various constructors.
		// Default.
		qbVector2();
		
		// With input data (std::vector).
		qbVector2(const std::vector<T> &inputData);
		// With input data (qbVector).
		qbVector2(const qbVector<T> &inputData);
		// With input data (qbVector2).
		qbVector2(const qbVector2<T> &inputData);
		// With input data as two separate values.
		qbVector2(const T x, const T y);
		
		// And the destructor.
		~qbVector2();	
		
		// Keep the GetNumDims() function for backwards compatibility.
		int GetNumDims() const;
		T GetElement(int index) const;
		void SetElement(int index, T value);		
		
		// Functions to perform computations on the vector.
		// Return the length of the vector.
		T norm();
		
		// Return a normalized copy of the vector.
		qbVector2<T> Normalized();
		
		// Normalize the vector in place.
		void Normalize();
		
		// Overloaded operators.
		qbVector2<T> operator+ (const qbVector2<T> &rhs) const;		
		qbVector2<T> operator- (const qbVector2<T> &rhs) const;
		qbVector2<T> operator* (const T &rhs) const;
		
		// Overload the assignment operator.
		qbVector2<T> operator= (const qbVector<T> &rhs);
		qbVector2<T> operator= (const std::vector<T> &rhs);
		qbVector2<T> operator= (const qbVector2<T> &rhs);
		
		// Friend functions.
		template <class U> friend qbVector2<U> operator* (const U &lhs, const qbVector2<U> &rhs);
		
		// Static functions.
		static T dot(const qbVector2<T> &a, const qbVector2<T> &b);
		
	public:
		T m_x;
		T m_y;
		
};

/* **************************************************************************************************
CONSTRUCTOR / DESTRUCTOR FUNCTIONS
/* *************************************************************************************************/
// The default constructor.
template <class T>
qbVector2<T>::qbVector2()
{
	m_x = static_cast<T>(0.0);
	m_y = static_cast<T>(0.0);
}

template <class T>
qbVector2<T>::qbVector2(const std::vector<T> &inputData)
{
	if (inputData.size() != 2)
		throw std::invalid_argument("Cannot assign std::vector to qbVector2 - assignment dimension mismatch.");
		
	m_x = inputData.at(0);
	m_y = inputData.at(1);
}

template <class T>
qbVector2<T>::qbVector2(const qbVector<T> &inputData)
{
	if (inputData.GetNumDims() != 2)
		throw std::invalid_argument("Cannot assign qbVector to qbVector2 - assignment dimension mismatch.");
		
	m_x = inputData.GetElement(0);
	m_y = inputData.GetElement(1);
}

template <class T>
qbVector2<T>::qbVector2(const qbVector2<T> &inputData)
{
	m_x = inputData.m_x;
	m_y = inputData.m_y;
}

template <class T>
qbVector2<T>::qbVector2(const T x, const T y)
{
	m_x = x;
	m_y = y;
}

template <class T>
qbVector2<T>::~qbVector2()
{
	// For now, we don't need to do anything in the destructor.
}

/* **************************************************************************************************
FUNCTIONS TO PERFORM COMPUTATIONS ON THE VECTOR
/* *************************************************************************************************/
// Compute the length of the vector, known as the 'norm'.
template <class T>
T qbVector2<T>::norm()
{		
	return sqrt((m_x*m_x) + (m_y*m_y));
}

// Return a normalized copy of the vector.
template <class T>
qbVector2<T> qbVector2<T>::Normalized()
{
	// Compute the vector norm.
	T vecNorm = this->norm();
	
	// Compute the normalized version of the vector.
	//qbVector<T> result(m_vectorData);
	qbVector2<T> result;
	result.m_x = m_x / vecNorm;
	result.m_y = m_y / vecNorm;

	return result;
}

// Normalize the vector in place.
template <class T>
void qbVector2<T>::Normalize()
{
	// Compute the vector norm.
	T vecNorm = this->norm();
	T denominator = static_cast<T>(1.0) / vecNorm;
	
	m_x = m_x / denominator;
	m_y = m_y / denominator;
}

/* **************************************************************************************************
OVERLOADED OPERATORS
/* *************************************************************************************************/
template <class T>
qbVector2<T> qbVector2<T>::operator+ (const qbVector2<T> &rhs) const
{
	qbVector2<T> result;
	result.m_x = m_x + rhs.m_x;
	result.m_y = m_y + rhs.m_y;
	
	return result;
}

template <class T>
qbVector2<T> qbVector2<T>::operator- (const qbVector2<T> &rhs) const
{
	qbVector2<T> result;
	result.m_x = m_x - rhs.m_x;
	result.m_y = m_y - rhs.m_y;
	
	return result;
}

template <class T>
qbVector2<T> qbVector2<T>::operator* (const T &rhs) const
{
	qbVector2<T> result;
	result.m_x = m_x * rhs;
	result.m_y = m_y * rhs;
	
	return result;
}

/* **************************************************************************************************
THE ASSIGNMENT (=) OPERATOR
/* *************************************************************************************************/
template <class T>
qbVector2<T> qbVector2<T>::operator= (const qbVector<T> &rhs)
{
	if (rhs.GetNumDims() != 2)
		throw std::invalid_argument("Cannot assign qbVector to qbVector2 - assignment dimension mismatch.");
	
	m_x = rhs.GetElement(0);
	m_y = rhs.GetElement(1);
	
	return *this;
}

template <class T>
qbVector2<T> qbVector2<T>::operator= (const std::vector<T> &rhs)
{
	if (rhs.size() != 2)
		throw std::invalid_argument("Cannot assign std::vector to qbVector2 - assignment dimension mismatch.");
	
	m_x = rhs.at(0);
	m_y = rhs.at(1);
	
	return *this;
}

template <class T>
qbVector2<T> qbVector2<T>::operator= (const qbVector2<T> &rhs)
{
	m_x = rhs.m_x;
	m_y = rhs.m_y;
	
	return *this;
}

/* **************************************************************************************************
FRIEND FUNCTIONS
/* *************************************************************************************************/
template <class T>
qbVector2<T> operator* (const T &lhs, const qbVector2<T> &rhs)
{
	// Perform scalar multiplication.
	qbVector2<T> result;
	result.m_x = lhs * rhs.m_x;
	result.m_y = lhs * rhs.m_y;
		
	return result;
}

/* **************************************************************************************************
STATIC FUNCTIONS
/* *************************************************************************************************/
template <class T>
T qbVector2<T>::dot(const qbVector2<T> &a, const qbVector2<T> &b)
{
	T cumulativeSum = (a.m_x * b.m_x) + (a.m_y * b.m_y);
	
	return cumulativeSum;
}

/* **************************************************************************************************
FUNCTIONS FOR COMPATIBILITY
/* *************************************************************************************************/
template <class T>
int qbVector2<T>::GetNumDims() const
{
	return 2;
}

template <class T>
T qbVector2<T>::GetElement(int index) const
{
	if (index == 0)
	{
		return m_x;
	}
	else if (index == 1)
	{
		return m_y;
	}
	else
	{
		throw std::invalid_argument("Attempt to get invalid element index.");
	}
	return 0;
}

template <class T>
void qbVector2<T>::SetElement(int index, T value)
{
	if (index == 0)
	{
		m_x = value;
	}
	else if (index == 1)
	{
		m_y = value;
	}
	else
	{
		throw std::invalid_argument("Attempt to set invalid element index.");
	}
}

#endif
