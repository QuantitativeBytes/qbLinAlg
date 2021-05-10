/* *************************************************************************************************

	TestCode_qbPCA
	
	  Code to test the principal component analysis functions.
	
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
#include "../qbPCA.h"

using namespace std;

int main()
{
	cout << "**********************************************" << endl;
	cout << "Testing principal component analysis code." << endl;
	cout << "**********************************************" << endl;
	cout << endl;
	
	{
		cout << "Testing with 100 observations of 3 variables:" << endl;
		
		// Read data from the .CSV file.
		string rowData;
		string number;
		stringstream rowDataStream;
		std::vector<double> testData;
		int numRows = 0;
		int numCols = 0;
		ifstream inputFile("qbPCATestData.csv");
		
		/* If the file open successfully then do stuff. */
		if (inputFile.is_open())
		{		
			cout << "Opened file successfully..." << endl;
			
			while (!inputFile.eof())
			{
			
				// Read the next line.
				getline(inputFile, rowData);
				
				// Loop through and extract the individual numbers.
				rowDataStream.clear();
				rowDataStream.str(rowData);
				
				if (numRows < 1)
					numCols = 0;
				
				while (rowDataStream.good())
				{
					getline(rowDataStream, number, ',');
					testData.push_back(atof(number.c_str()));
					
					if (numRows < 1)
						numCols++;
				}
				numRows++;
			}
						
			// Close the file.
			inputFile.close();
			
			// For some reason, the above reads one line too many.
			numRows--;
			testData.pop_back();
					
			cout << "Completed reading file..." << endl;
			cout << "Read " << numRows << " observations of " << numCols << " variables." << endl;
			cout << "Constituting " << testData.size() << " elements in total." << endl;
			
			// Form into a matrix.
			qbMatrix2<double> X (numRows, numCols, testData);
			
			// Compute the covariance matrix.
			std::vector<double> columnMeans = qbPCA::ComputeColumnMeans(X);
			qbMatrix2<double> X2 = X;
			qbPCA::SubtractColumnMeans(X2, columnMeans);
			
			qbMatrix2<double> covX = qbPCA::ComputeCovariance(X2);
			cout << endl;
			cout << "Giving the covariance matrix as: " << endl;
			covX.PrintMatrix();
			
			// Compute the eigenvectors.
			qbMatrix2<double> eigenvectors;
			int testResult = qbPCA::ComputeEigenvectors(covX, eigenvectors);
			cout << endl;
			cout << "And the eigenvectors as: " << endl;
			eigenvectors.PrintMatrix();
			
			// Test the overall function.
			cout << endl;
			cout << "Testing overall function..." << endl;
			qbMatrix2<double> eigenvectors2;
			int testResult2 = qbPCA::qbPCA(X, eigenvectors2);
			cout << "testResult2 = " << testResult2 << endl;
			cout << "And the final eigenvectors are:" << endl;
			eigenvectors2.PrintMatrix();
			
			// Test dimensionality reduction.
			cout << endl;
			cout << "Testing dimensionality reduction." << endl;
			cout << "Starting with X which has " << X.GetNumRows() << " rows and " << X.GetNumCols() << " columns." << endl;
			cout << endl;
			cout << "Using only the first two principal components:" << endl;
			qbMatrix2<double> V, part2;
			eigenvectors.Separate(V, part2, 2);
			V.PrintMatrix(8);
			cout << endl;
			
			qbMatrix2<double> newX = (V.Transpose() * X.Transpose()).Transpose();
			cout << "Result has " << newX.GetNumRows() << " rows and " << newX.GetNumCols() << " columns." << endl;
			
			// Open a file for writing
			ofstream outputFile("qbPCATestData_Reduced.csv");
			if (outputFile.is_open())
			{
				for (int i=0; i<newX.GetNumRows(); ++i)
				{
					outputFile << newX.GetElement(i, 0) << "," << newX.GetElement(i, 1) << endl;
				}
				outputFile.close();
			}			
			
		}
	}
}
