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
  cout << std::fixed << std::setprecision(3) << inputVector.GetElement(row) << endl;
	}    
}

int main()
{
	cout << "**********************************************" << endl;
	cout << "Testing eigenvalue and eigenvector code." << endl;
	cout << "**********************************************" << endl;
	cout << endl;
	
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
	qbEIG_PIt(testMatrix, eigenValue, eigenVector);
	
	cout << "Eigenvector: " << endl;
	PrintVector(eigenVector);
	cout << "Eigenvalue = " << eigenValue << "." << endl;
	cout << endl;
	
	return 0;
}
