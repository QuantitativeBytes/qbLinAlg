/* *************************************************************************************************

	TestCode_qbQR
	
	  Code to test the QR decomposition code.
	
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
#include "../qbQR.h"

using namespace std;

int main()
{
	cout << "**********************************************" << endl;
	cout << "Testing QR decomposition code." << endl;
	cout << "**********************************************" << endl;
	cout << endl;

	{	
		cout << "Testing with simple 3x3 matrix:" << endl;
	
		std::vector<double> simpleData = {0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
		qbMatrix2<double> testMatrix(3, 3, simpleData);
		
		testMatrix.PrintMatrix();
		
		cout << endl;
		cout << "Computing QR decomposition..." << endl;
		qbMatrix2<double> Q (3,3);
		qbMatrix2<double> R (3,3);
		int status = qbQR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.PrintMatrix();
		cout << endl;
		cout << "Q = " << endl;
		Q.PrintMatrix();
		cout << endl;
	}
	
	{
		cout << "Testing with simple 4x4 matrix:" << endl;
		std::vector<double> simpleData = {1.0, 5.0, 3.0, 4.0, 7.0, 8.0, 2.0, 9.0, 7.0, 3.0, 2.0, 1.0, 9.0, 3.0, 5.0, 7.0};
		qbMatrix2<double> testMatrix(4, 4, simpleData);
		testMatrix.PrintMatrix();
		cout << endl;
		cout << "Computing QR decomposition..." << endl;
		qbMatrix2<double> Q (4,4);
		qbMatrix2<double> R (4,4);
		int status = qbQR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.PrintMatrix();
		cout << endl;
		cout << "Q = " << endl;
		Q.PrintMatrix();
		cout << endl;		
	}
	
	return 0;
}
