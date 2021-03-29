/* *************************************************************************************************

	TestCode_qbEIG
	
	  Code to test the eigenvector and eigenvalue code.
	
	*** INPUTS ***
	
	None
															
	*** OUTPUTS ***
	
	INT				Flag indicating success or failure of the process.

	Created as part of the qbLinAlg linear algebra library, which is intended to be primarily for
	educational purposes. For more details, see the corresponding videos on the QuantitativeBytes
	YouTube channel at:
	
	www.youtube.com/c/QuantitativeBytes								

	************************************************************************************************* */

#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <fstream>

#include "../qbMatrix.h"
#include "../qbVector.h"
#include "../qbEIG.h"

using namespace std;

// A simple function to print a matrix to stdout.
template <class T>
void PrintMatrix(qbMatrix2<T> matrix)
{
	int nRows = matrix.GetNumRows();
	int nCols = matrix.GetNumCols();
	for (int row = 0; row<nRows; ++row)
  {
	  for (int col = 0; col<nCols; ++col)
    {
	    cout << std::fixed << std::setprecision(3) << matrix.GetElement(row, col) << "  ";
    }
	cout << endl;
	}    
}

// A simple function to print a vector to stdout.
template <class T>
void PrintVector(qbVector<T> inputVector)
{
	int nRows = inputVector.GetNumDims();
	for (int row = 0; row<nRows; ++row)
  {
  cout << std::fixed << std::setprecision(6) << inputVector.GetElement(row) << endl;
	}    
}

int main()
{
	cout << "**********************************************" << endl;
	cout << "Testing eigenvalue and eigenvector code." << endl;
	cout << "Power Iteration Method." << endl;
	cout << "**********************************************" << endl;
	cout << endl;
	
	{
		cout << "Testing with simple 3x3 matrix:" << endl;
		
		//std::vector<double> simpleData = {1.0, 2.0, 3.0, 4.0};
		//qbMatrix2<double> testMatrix(2, 2, simpleData);
	
		std::vector<double> simpleData = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
		qbMatrix2<double> testMatrix(3, 3, simpleData);
		
		PrintMatrix(testMatrix);
		
		cout << endl;
		cout << "Computing eigenvector and eigenvalue..." << endl;
		double eigenValue;
		qbVector<double> eigenVector;
		qbEIG_PIt<double>(testMatrix, eigenValue, eigenVector);
		
		cout << "Eigenvector: " << endl;
		PrintVector(eigenVector);
		cout << "Eigenvalue = " << eigenValue << "." << endl;
		cout << endl;
	}
	
	cout << "**********************************************" << endl;
	cout << "Testing eigenvalue and eigenvector code." << endl;
	cout << "Inverse-Power Iteration Method." << endl;
	cout << "**********************************************" << endl;
	cout << endl;	
	
	{
		cout << "Testing with a simple 3x3 matrix:" << endl;
		std::vector<double> simpleData = {0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
		qbMatrix2<double> testMatrix(3, 3, simpleData);
		PrintMatrix(testMatrix);
		cout << endl;
		
		// Define the eigenvalues.
		std::vector<double> eigenValues = {1.5962551, -0.37253087, 0.02627577};
		
		// Setup a vector for the eigenvector.
		qbVector<double> eigenVector(3);
		
		// Loop through each eigenvalue.
		for (auto currentValue : eigenValues)
		{
			cout << "Estimated eigenvector for eigenvalue " << currentValue << " = " << endl;
			int returnStatus = qbInvPIt<double>(testMatrix, currentValue, eigenVector);
			PrintVector(eigenVector);
			
			if (returnStatus == QBEIG_MAXITERATIONSEXCEEDED)
				cout << "*** Maximum iterations exceeded ***" << endl;
			
			cout << endl;
		}
		
	}
	
	cout << "**********************************************" << endl;
	cout << "Testing eigenvalue and eigenvector code." << endl;
	cout << "New IsSymmetric() function in qbMatrix2 class." << endl;
	cout << "**********************************************" << endl;
	cout << endl;			
	
	{
		cout << "Testing with a simple symmetric matrix." << endl;
		std::vector<double> simpleData = {4.0, -7.0, 6.0, -7.0, -2.0, 13.0, 6.0, 13.0, 5.0};
		qbMatrix2<double> testMatrix(3, 3, simpleData);
		PrintMatrix(testMatrix);
		cout << endl;
		
		cout << "Matrix is symmetric: ";
		if (testMatrix.IsSymmetric())
			cout << "True." << endl;
		else
			cout << "False." << endl;
	}
	
	{
		cout << "Testing with a simple non-symmetric matrix." << endl;
		std::vector<double> simpleData = {4.0, -7.0, 6.0, 7.0, -2.0, 13.0, 6.0, 13.0, 5.0};
		qbMatrix2<double> testMatrix(3, 3, simpleData);
		PrintMatrix(testMatrix);
		cout << endl;
		
		cout << "Matrix is symmetric: ";
		if (testMatrix.IsSymmetric())
			cout << "True." << endl;
		else
			cout << "False." << endl;
	}		
	
	cout << endl << endl;	
	cout << "**********************************************" << endl;
	cout << "Testing eigenvalue and eigenvector code." << endl;
	cout << "Eigenvalues by QR decomposition." << endl;
	cout << "**********************************************" << endl;
	cout << endl;		
	
	{
		cout << "Testing with a simple (non-symmetric) 3x3 matrix:" << endl;
		std::vector<double> simpleData = {0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
		qbMatrix2<double> testMatrix(3, 3, simpleData);
		PrintMatrix(testMatrix);
		cout << endl;
		
		// Compute the eigenvalues.
		std::vector<double> eigenValues;
		int returnStatus = qbEigQR(testMatrix, eigenValues);
		
		if (returnStatus == QBEIG_MAXITERATIONSEXCEEDED)
			cout << ">>> Maximum iterations exceeded <<<" << endl;
			
		if (returnStatus == QBEIG_MATRIXNOTSYMMETRIC)
			cout << ">>> Matrix not symmetric. <<<" << endl;	
		
		// Display the eigenvalues.
		cout << "The estimated eigenvalues are:" << endl;
		for (auto currentValue : eigenValues)
			cout << std::setprecision(6) << currentValue << " ";
			
		cout << endl << endl;	
	}
	
	{
		cout << "Testing with a simple (symmetric) 3x3 matrix:" << endl;
		std::vector<double> simpleData = {6.0, 5.5, -1.0, 5.5, 1.0, -2.0, -1.0, -2.0, -3.0};
		qbMatrix2<double> testMatrix(3, 3, simpleData);
		PrintMatrix(testMatrix);
		cout << endl;
		
		// Compute the eigenvalues.
		std::vector<double> eigenValues;
		int returnStatus = qbEigQR(testMatrix, eigenValues);
		
		if (returnStatus == QBEIG_MAXITERATIONSEXCEEDED)
			cout << ">>> Maximum iterations exceeded <<<" << endl;
			
		if (returnStatus == QBEIG_MATRIXNOTSYMMETRIC)
			cout << ">>> Matrix not symmetric. <<<" << endl;	
		
		// Display the eigenvalues.
		cout << "The estimated eigenvalues are:" << endl;
		for (auto currentValue : eigenValues)
			cout << std::setprecision(6) << currentValue << " ";
			
		cout << endl << endl;
			
		// Setup a vector for the eigenvector.
		qbVector<double> eigenVector(3);
		
		// Loop through each eigenvalue.
		for (auto currentValue : eigenValues)
		{
			cout << "Estimated eigenvector for eigenvalue " << currentValue << " = " << endl;
			int returnStatus = qbInvPIt<double>(testMatrix, currentValue, eigenVector);
			PrintVector(eigenVector);
			
			if (returnStatus == QBEIG_MAXITERATIONSEXCEEDED)
				cout << "*** Maximum iterations exceeded ***" << endl;
			
			cout << endl;
		}			
			
		cout << endl << endl;	
	}	
	
	{
		cout << "Testing with an example that should have complex eigenvalues:" << endl;
		std::vector<double> simpleData = {4.0, -6.0, 8.0, 7.0, 9.0, -5.0, 9.0, -6.0, -4.0};
		qbMatrix2<double> testMatrix(3, 3, simpleData);
		PrintMatrix(testMatrix);
		cout << endl;
		
		// Compute the eigenvalues.
		std::vector<double> eigenValues;
		int returnStatus = qbEigQR(testMatrix, eigenValues);
		
		if (returnStatus == QBEIG_MATRIXNOTSYMMETRIC)
			cout << ">>> Matrix not symmetric. <<<" << endl;		
		
		if (returnStatus == QBEIG_MAXITERATIONSEXCEEDED)
			cout << ">>> Maximum iterations exceeded <<<" << endl;			
		
		// Display the eigenvalues.
		cout << "The estimated eigenvalues are:" << endl;
		for (auto currentValue : eigenValues)
			cout << std::setprecision(6) << currentValue << " ";
			
		cout << endl << endl;	
		
		// Setup a vector for the eigenvector.
		qbVector<double> eigenVector(3);		
		
		// Loop through each eigenvalue.
		for (auto currentValue : eigenValues)
		{
			cout << "Estimated eigenvector for eigenvalue " << currentValue << " = " << endl;
			int returnStatus = qbInvPIt<double>(testMatrix, currentValue, eigenVector);
			PrintVector(eigenVector);
			
			if (returnStatus == QBEIG_MAXITERATIONSEXCEEDED)
				cout << "*** Maximum iterations exceeded ***" << endl;
			
			cout << endl;
		}
		
		cout << endl << endl;
	}
	
	return 0;
}
