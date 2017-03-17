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

/*TODO 1: update functions to use arrays.
  TODO 2: validate inputs: http://www.cplusplus.com/forum/beginner/13044/#msg62827
          Can either use template gist or notebook loop tests.
  TODO 3: 


 */
 
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>

using namespace std; // TODO: remove this and update to use std::


struct equation_t {
  int coeffsCount = 3;  // Override Counts for other degree polynomials
  int rootsCount = 2;  
  double* coeffs;  // pointer to coefficient array
  double* roots;  // pointer to roots array
  bool flag;  // if true then real roots exist
};

// Creates arrays for coeffs and roots.
void equation_init(equation_t&);

// Reads and returns valid coefficients from stdin.
void readCoeffs(equations_t&); 

// Calculates and returns the discriminant.
double discr(equations_t&);

// If the discriminant is zero or greater, the roots are computed
// and stored in global variables, then returns true.
bool equSolver(equations_t&);

// Appends real roots to the file. Otherwise prints a message
// to stdout that no real roots exists.
void outResults(equations_t&, ofstream&);

// equation_cleanup frees memory resourses before program exits
void equation_cleanup(equations_t&, int);



int main(int argc, char* argv[]) {
  // local constants
  constexpr char* OUTPUT_FILE = "results.dat";

  // local variables
  int number = (argc > 1)?atoi(argv[1]):0;  // number of quadratic formulas to calculate
  
  // struct array to hold all the equations calculated
  equations_t quadratic[number];
  
  ofstream outStream;
  outStream.open(OUTPUT_FILE);  // if file exists, it is overwritten
  
  // Each iteration calculates one quadratic equation.
  for (int i=0; i < number; i++){
    // initialize this quadratic equation
    equation_init(quadratic[i]);
    
    // Operator enters values for the coefficients.
    readCoeffs(quadratic[i]);
    
    // Determines if there are real roots. If so calculates the roots. 
    equSolver(quadratic[i]);
    
    // Directs output to either the file or stdout respectively, based on flag.
    outResults(quadratic[i], outStream);
  }
  outStream.close
  equation_cleanup(quadratic, number)
  return 0;
}

// Assigns address of new arrays to coeffs and roots pointers
void equation_init(equation_t& eq){
  // These arrays are allocated from the heap and exist until deleted.
  // Next version I plan to handle that with a class destructor.
  eq.coeffs = new double[eq.coeffsCount];  
  eq.roots = new double[eq.rootsCount];
  return;
}

// Operator inputs coefficients.  If zero is entered for coeffiecient a, an
// error message is displayed requesting a new entry.
void readCoeffs(double coeffs, int coeffsSize){

  while (true){  // Runs ad-infinitum until break condition is met.
    std::std::cout << "\nEnter coefficient a: "; 
    std::cin >> *coeffs;
    if (*coeffs) break;  // a must not equal zero 
    else {  // operator entered zero for coefficient a
      std::cout << "\nInvalid entry. Please enter a non-zero "\
           "value for a." << std::endl;
    }
  }
  std::cout << "\nEnter coefficient b: "; 
  std::cin >> *++coeffs;
  std::cout << "\nEnter coefficient c: "; 
  std::cin >> *++coeffs;
  return; 
}

// Computes and returns the discriminant.
double discr(double coeffs[], int coeffsSize){
  return (coeffs[1]*coeffs[1]-4*coeffs[0]*coeffs[2]);
}

// Gets the discriminant and if it's greater than or equal to 0 
// computes the roots and returns true.
bool equSolver(double coeffs[], int coeffsSize, double roots[], int rootsSize){
  double compDisc = discr(coeffs, coeffsSize);
  
  if (compDisc >= 0){
    roots[0] = (-coeffs[1] + sqrt(compDisc))/(2*coeffs[0]);
    roots[1] = (-coeffs[1] - sqrt(compDisc))/(2*coeffs[0]);
  }
  
  return (compDisc >= 0)?true:false; // If roots exists true is returned.
}

void outResults(double coeffs[], int coeffsSize, bool rootsExist, ofstream& outStream){
  if (rootsExist){
    // Results are appended to the opened file.
    outStream << "Quadratic equation with the following coefficients:";
    outStream << std::endl << "a: " << a << "; b: " << b << "; c: " << c << std::endl;
    outStream << "has the following roots" << std::endl << "Root1: " << root1 ;
    outStream << "; Root2: " << root2 << ";" << std::endl << std::endl;
  } else {
    std::cout << "Quadratic equation with the following coefficients:" << std::endl;
    std::cout << "a: " << a << "; b: " << b << "; c: " << c << std::endl;
    std::cout << "has no roots in the real domain." << std::endl << std::endl;
  }
  return;
}

void equation_cleanup(equations_t& eq, int number){
  for (int i=0; i<number;i++){
    
    // I don't know if I have to dereference these for this to work
    delete eq[i].coeffs; 
    delete eq[i].roots;
    
    // not really needed but good practice i think
    eq[i].coeffs = nullptr;
    eq[i].roots = nullptr;
  }
}