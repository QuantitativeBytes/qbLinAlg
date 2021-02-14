/* *************************************************************************************************

	TestCode_qbLinearSolve
	
	 Code to test solving systems of linear equations using the qbLinearSolve function. Also tests
	computation of the matrix rank in situations where Gaussian eliminiation is possible and also
	where it is not.
	
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

#include "../qbMatrix.h"
#include "../qbVector.h"
#include "../qbLinSolve.h"

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
	cout << "Code to test qbMatrix2" << endl;
	cout << "Testing conversion of matrix to row echelon form." << endl;
	cout << endl;
	
	// Generate a matrix to test things with.
  std::vector<double> simpleData = {1.0, 3.0, -1.0, 13.0, 4.0, -1.0, 1.0, 9.0, 2.0, 4.0, 3.0, -6.0};
  qbMatrix2<double> testMatrix(3, 4, simpleData);
  
  cout << "Original matrix:" << endl;
  PrintMatrix(testMatrix);
  cout << endl;
  
  // Convert to row echelon form.
  qbMatrix2<double> rowEchelonMatrix = testMatrix.RowEchelon();
  
  cout << "Converted to row echelon form:" << endl;
  PrintMatrix(rowEchelonMatrix);
  cout << endl;
  
  // Define another matrix as the first part of our system of linear equations.
  std::vector<double> simpleData2 = {1.0, 3.0, -1.0, 4.0, -1.0, 1.0, 2.0, 4.0, 3.0};
  qbMatrix2<double> aMat(3, 3, simpleData2);
  cout << "We setup the equations in the form of Ax = b, where A = " << endl;
  PrintMatrix(aMat);
  cout << endl;
  
  // Define a vector to hold the RHS of our system of linear equations.
  std::vector<double> vectorData {13.0, 9.0, -6.0};
  qbVector<double> bVec {vectorData};
  cout << "And b = " << endl;
  PrintVector(bVec);
  cout << endl;
  
  // Call the qbLinSolve function.
  qbVector<double> testResult(3);
  int test = qbLinSolve<double>(aMat, bVec, testResult);
  cout << "And the final result is:" << endl;
  PrintVector(testResult);
  cout << endl;
  
  // ***************************************************************************************************
  // Try some random tests.
  std::random_device myRandomDevice;
  std::mt19937 myRandomGenerator(myRandomDevice());
  std::uniform_real_distribution<double> myDistribution(-25.0,25.0);
  
  int numUnknowns = 10;
  
  std::vector<double> coefficientData;
  std::vector<double> unknownData;
  // Populate the coefficient data.
  for (int i=0; i<(numUnknowns * numUnknowns); ++i)
  {
  	double randomNumber = myDistribution(myRandomGenerator);
  	coefficientData.push_back(randomNumber);
  }
  cout << "A random coefficient matrix = " << endl;
  qbMatrix2<double> coefficientMatrix(numUnknowns, numUnknowns, coefficientData);
  PrintMatrix(coefficientMatrix);
  cout << endl;
  
  cout << "And the random unknown values = " << endl;
  for (int i=0; i<numUnknowns; ++i)
  {
  	double randomNumber = myDistribution(myRandomGenerator);
  	unknownData.push_back(randomNumber);
  }
  qbVector<double> unknownVector {unknownData};
  PrintVector(unknownVector);
  cout << endl;
  
  cout << "Compute the equation results = " << endl;
  qbVector<double> systemResult = coefficientMatrix * unknownVector;
  PrintVector(systemResult);
  cout << endl;
  
  cout << "Attempt to solve the linear system..." << endl;
  qbVector<double> compSolution(numUnknowns);
  int compTest = qbLinSolve<double>(coefficientMatrix, systemResult, compSolution);
  PrintVector(compSolution);
  cout << endl;
  
  cout << "And compare the actual result with the computed solution..." << endl;
  qbVector<double> errorVector = unknownVector - compSolution;
  PrintVector(errorVector);
  cout << endl;
  
  // ***************************************************************************************************
  // Test computation of the matrix rank.
  cout << "***************************************************************" << endl;
  cout << "Testing computation of matrix rank" << endl;
  cout << "***************************************************************" << endl;
  cout << endl;
  cout << "Testing with a solvable system:" << endl;
  PrintMatrix(aMat);
  cout << "Rank = " << aMat.Rank() << endl;
  cout << endl;
  cout << "Row echelon form:" << endl;
  qbMatrix2<double> aMatRowEchelon = aMat.RowEchelon();
  PrintMatrix(aMatRowEchelon);
  cout << endl;
  
  // Test the condition when Gaussian elmination fails.
  cout << "***************************************************************" << endl;
  cout << "Testing the condition when Gaussian elimination fails" << endl;
  cout << "***************************************************************" << endl;
  std::vector<double> geFailData = {0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0};
  qbMatrix2<double> geFailMatrix (3, 3, geFailData);
  cout << "Testing with:" << endl;
  PrintMatrix(geFailMatrix);
  cout << endl;
  cout << "Attempt to perform Gaussian elimination on this gives:" << endl;
  qbMatrix2<double> geFailResult = geFailMatrix.RowEchelon();
  PrintMatrix(geFailResult);
  cout << endl;
  cout << "Attempt to compute the rank gives:" << endl;
  cout << "Rank = "<< geFailMatrix.Rank() << endl;
  cout << endl;
  cout << "Testing with a zero matrix:" << endl;
  std::vector<double> geFailData2 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  qbMatrix2<double> geFailMatrix2 (3, 3, geFailData2);
  PrintMatrix(geFailMatrix2);
  cout << endl;
  cout << "Rank = " << geFailMatrix2.Rank() << endl;
  cout << endl;
  cout << "Testing with a larger example:" << endl;
  PrintMatrix(coefficientMatrix);
  cout << endl;
  cout << "The rank is " << coefficientMatrix.Rank() << endl;
  cout << endl;
  
  std::vector<double> geFailData3 = {1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0, 2.0, 4.0, 6.0, 8.0, 8.0, 6.0, 4.0, 2.0};
  cout << "Testing with 4x4:" << endl;
  cout << endl;
  qbMatrix2<double> geFailMatrix3 (4, 4, geFailData3);
  PrintMatrix(geFailMatrix3);
  cout << "Rank = " << geFailMatrix3.Rank() << endl;
  cout << endl;
  
  // Test the two possible conditions with no solution
  cout << "***************************************************************" << endl;
  cout << "Testing the two possible conditions with no solution" << endl;
  cout << "***************************************************************" << endl;  
  cout << endl;
  {
  	// Setup a system with an infinite number of solutions.
  	cout << "A system with an infinite number of solutions:" << endl;
  	std::vector<double> aMatData = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0};
  	std::vector<double> bVecData = {0.0, 0.0, 0.0};
  	qbMatrix2<double> aMat (3, 3, aMatData);
  	qbVector<double> bVec {bVecData};
  	qbVector<double> solution(3);
  	int test = qbLinSolve<double>(aMat, bVec, solution);
  	if (test > 0)
	  	PrintVector<double>(solution);
	  else
	  	cout << "Error condition: " << test << endl;
  }
  cout << endl;
  {
  	// Setup a system with no solutions.
  	cout << "A system with no solutions:" << endl;
  	std::vector<double> aMatData = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
  	std::vector<double> bVecData = {0.0, -1.0, 1.0};
  	qbMatrix2<double> aMat (3, 3, aMatData);
  	qbVector<double> bVec {bVecData};
  	qbVector<double> solution(3);
  	int test = qbLinSolve<double>(aMat, bVec, solution);
  	if (test > 0)
	  	PrintVector<double>(solution);
	  else
	  	cout << "Error condition: " << test << endl;
  }  
  
	return 0;
}   
