/***************************************************************************
 file   : boatwright3.cpp
 title  : COP 2001 Assignment #3 - 2017
 author : Ivan Boatwright
 email  : ijboatwright5153@eagle.fgcu.edu
 version: 3.0 2/28/17
 
   Solves quadratic equations of the form:
        a*x^2 + b*x + c = 0
   The number of equations to calculate is passed from the commandline.
   The values for a, b, and c are entered by the operator when prompted.
   If zero is entered for the 'a' variable an error message is displayed.
   If roots exist the results are printed to a file. If no roots exist a
   message is printed to stdout. 
 ***************************************************************************/

//TODO: update prototypes and functions to use arrays
 
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>

using namespace std;

// Reads and returns valid coefficients from stdin.
void readCoeffs(double&, double&, double&); 

// Calculates and returns the discriminant.
double discr(double, double, double);

// If the discriminant is zero or greater, the roots are computed
// and stored in global variables, then returns true.
bool equSolver(double, double, double);

// Appends real roots to the file. Otherwise prints a message
// to stdout that no real roots exists.
void outResults(double, double, double, bool, ofstream&);

int main(int argc, char* argv[]) {
  // local constants
  const char* OUTPUT_FILE = "results.dat";

  // local variables
  double coeffs[3];  // coefficients
  double roots[2];  // roots
  bool flag;  // if true then real roots exist
  int number = (argc > 1)?atoi(argv[1]):0;  // number of quadratic formulas to calculate
  
  ofstream outStream;
  outStream.open(OUTPUT_FILE);  // if file exists, it is overwritten
  
  // Each iteration calculates one quadratic equation.
  for (int i=0; i < number; i++){
    // Operator enters values for the coefficients.
    readCoeffs(a, b, c);
    
    // Determines if there are real roots. If so calculates the roots. 
    flag = equSolver(a, b, c);
    
    // Directs output to either the file or stdout respectively, based on flag.
    outResults(a, b, c, flag, outStream);
  }
  outStream.close();
  return 0;
}

// Operator inputs coefficients.  If zero is entered for coeffiecient a, an
// error message is displayed requesting a new entry.
void readCoeffs(double& a, double& b, double& c){
  while (true){  // Runs ad-infinitum until break condition is met.
    cout << "\nEnter coefficient a: "; 
    cin >> a;
    if (a) break;  // a must not equal zero 
    else {  // operator entered zero for coefficient a
      cout << "\nInvalid entry. Please enter a non-zero "\
           "value for a." << endl;
    }
  }
  cout << "\nEnter coefficient b: "; 
  cin >> b;
  cout << "\nEnter coefficient c: "; 
  cin >> c;
  return; 
}

// Computes and returns the discriminant.
double discr(double a, double b, double c){
  return (b*b-4*a*c);
}

// Gets the discriminant and if it's greater than or equal to 0 
// computes the roots and returns true.
bool equSolver(double a, double b, double c){
  double compDisc = discr(a, b, c);
  
  if (compDisc >= 0){
    root1 = (-b + sqrt(compDisc))/(2*a);
    root2 = (-b - sqrt(compDisc))/(2*a);
  }
  
  return (compDisc >= 0)?true:false; // If roots exists true is returned.
}

void outResults(double a, double b, double c, bool rootsExist, ofstream& outStream){
  if (rootsExist){
    // Results are appended to the opened file.
    outStream << "Quadratic equation with the following coefficients:";
    outStream << endl << "a: " << a << "; b: " << b << "; c: " << c << endl;
    outStream << "has the following roots" << endl << "Root1: " << root1 ;
    outStream << "; Root2: " << root2 << ";" << endl << endl;
  } else {
    cout << "Quadratic equation with the following coefficients:" << endl;
    cout << "a: " << a << "; b: " << b << "; c: " << c << endl;
    cout << "has no roots in the real domain." << endl << endl;
  }
  return;
}
