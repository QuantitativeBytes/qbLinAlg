/* *************************************************************************************************

	qbMatrix_DeterminantTest
	
	Code to test computation of the determinant using the qbMatrix2 class.
	
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
#include <fstream>

#include "../qbMatrix.h"

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

int main()
{

	cout << "Testing implementation of determinant calculation." << endl;
	cout << endl;
	
	cout << "Generate a test matrix." << endl;
	double testData[9] = {2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 3.0, 1.0};
	qbMatrix2<double> testMatrix(3, 3, testData);
	PrintMatrix(testMatrix);
	cout << endl;
	
	cout << "Extract sub-matrix for element (0,0)" << endl;
	qbMatrix2<double> minor1 = testMatrix.FindSubMatrix(0,0);
	PrintMatrix(minor1);
	cout << endl;
	
	cout << "Extract sub-matrix for element (0,1)" << endl;
	qbMatrix2<double> minor2 = testMatrix.FindSubMatrix(0,1);
	PrintMatrix(minor2);
	cout << endl;
	
	cout << "Extract sub-matrix for element (0,2)" << endl;
	qbMatrix2<double> minor3 = testMatrix.FindSubMatrix(0,2);
	PrintMatrix(minor3);
	cout << endl;
	
	cout << "Extract sub-matrix for element (1,1)" << endl;
	qbMatrix2<double> minor4 = testMatrix.FindSubMatrix(1,1);
	PrintMatrix(minor4);
	cout << endl;	

	cout << "Test with a larger matrix." << endl;
	double testData2[25] =
		{2.0, 3.0, 4.0, 5.0, 6.0,
		 1.0, 2.0, 3.0, 4.0, 5.0,
		 9.0, 5.0, 3.0, 2.0, 6.0,
		 2.0, 4.0, 6.0, 5.0, 1.0,
		 1.0, 7.0, 5.0, 2.0, 3.0};
	qbMatrix2<double> testMatrix2(5, 5, testData2);
	PrintMatrix(testMatrix2);
	cout << endl;

	cout << "Extract sub-matrix for element (0,0)" << endl;
	qbMatrix2<double> minor5 = testMatrix2.FindSubMatrix(0,0);
	PrintMatrix(minor5);
	cout << endl;
	
	cout << "Extract sub-matrix for element (0,1)" << endl;
	qbMatrix2<double> minor6 = testMatrix2.FindSubMatrix(0,1);
	PrintMatrix(minor6);
	cout << endl;
	
	cout << "Extract sub-matrix for element (0,2)" << endl;
	qbMatrix2<double> minor7 = testMatrix2.FindSubMatrix(0,2);
	PrintMatrix(minor7);
	cout << endl;
	
	cout << "Extract sub-matrix for element (1,1)" << endl;
	qbMatrix2<double> minor8 = testMatrix2.FindSubMatrix(1,1);
	PrintMatrix(minor8);
	cout << endl;	

	cout << "Test determinant of 3x3 matrix:" << endl;
	cout << testMatrix.Determinant() << endl;
	cout << endl;

	cout << "Test determinant of 5x5 matrix:" << endl;
	cout << testMatrix2.Determinant() << endl;
	cout << endl;
	
	cout << "Test determinant of a singular matrix:" << endl;
	double testData3[9] =
		{1.0, 1.0, 1.0,
		 0.0, 1.0, 0.0,
		 1.0, 0.0, 1.0};	
	qbMatrix2<double> testMatrix3(3, 3, testData3);
	PrintMatrix(testMatrix3);
	cout << endl;
	cout << "Determinant = " << testMatrix3.Determinant() << endl;
	cout << endl;

	return 0;
}
