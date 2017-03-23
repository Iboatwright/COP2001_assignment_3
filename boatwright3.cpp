/***************************************************************************
 file   : boatwright3.cpp
 title  : COP 2001 Assignment #3 - 2017
 author : Ivan Boatwright
 email  : ijboatwright5153@eagle.fgcu.edu
 version: 3.14 3/22/17
 
   Solves quadratic equations of the form:
        a*x^2 + b*x + c = 0
   The number of equations to calculate is passed from the commandline.
   Only a positive integer is accepted from the commandline. Otherwise the
     program terminates early.
   The values for a, b, and c are entered by the operator when prompted.
    If zero is entered for the 'a' variable an error message is displayed.
    If non-numeric values are entered for a, b, or c an error message is
      displayed.
   If roots exist the results are printed to a file. If no roots exist a
   message is printed to stdout. 
 ***************************************************************************/
 
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <limits>

// defaults to a quadratic 
struct equation_t {
  int coeffsCount = 3;  // Override Counts for other degree polynomials
  int rootsCount = 2;  
  double* coeffs;  // pointer to coefficient array
  double* roots;  // pointer to roots array
  bool rootsExist = false;  // if true then real roots exist
};

// commandline argument validation (must be a positive integer)
int validateNumber(int, char*, int);

// Creates new[] equation arrays for coeffs and roots. Resources should be
//  released before program exits.
void equation_init(equation_t&);

// Requests user to input the coefficients
void coeffInput(double&, char);

// Reads and returns valid coefficients from stdin.
void readCoeffs(equation_t&); 

// Calculates and returns the discriminant.
double discr(double*);

// If the discriminant is zero or greater, the roots are computed,
// stored in the equation object, and rootsExist is set to true.
void equSolver(equation_t&);

// Appends real roots to the file. Otherwise prints a message
// to stdout that no real roots exists.
void outResults(equation_t&, std::ofstream&);

// releases memory resourses before program exits
void resource_cleanup(equation_t*, int);

int main(int argc, char* argv[]) {
  // local constants
  const char* OUTPUT_FILE = "results.dat";
  const int max_number_size = 10000; // prevents stack overflow

  // number of equations to calculate
  int number = validateNumber(argc, argv[1], max_number_size);  

  // struct array to hold all the equations calculated
  equation_t quadratic[number];
  
  std::ofstream outStream;
  outStream.open(OUTPUT_FILE);  // if file exists, it is overwritten
  
  // Each iteration calculates one quadratic equation.
  for (int i=0; i < number; i++){
    // Initializes this quadratic equation instance.
    equation_init(quadratic[i]);
    
    // Operator enters values for the coefficients.
    readCoeffs(quadratic[i]);
    
    // Determines if there are real roots. If so calculates the roots. 
    equSolver(quadratic[i]);
    
    // Directs output to the file or stdout respectively, based on rootsExist.
    outResults(quadratic[i], outStream);
  }
  
  outStream.close();
  resource_cleanup(quadratic, number);
  
  // Friendly end of program message.
  if ( number > 1){
    std::cout << "( " << number << " ) equation" << ((number == 1)?" ":"s ")
                    << "calculated. Have a nice day!" << std::endl;
  } 
  
  return 0;
}

// For invalid input, a message is displayed and 0 is returned
int validateNumber(int argc, char* argv1, int max_number_size){
  // atoi returns 0 for anything it can't convert to an int
  int number = (argc > 1)?std::atoi(argv1):0;
  if (argc > 1 && number < 1) { // not a positive integer
    number = 0; // discards any negative values 
    std::cout << "Program terminated.\n" << argv1 << " is an invalid "\
              "number of equations.\nPlease enter a positive integer "\
              "for the number of equations to calculate\n\tand run the program again." << std::endl;
  } else if (argc == 1) { // nothing entered on the commandline
    std::cout << "Program terminated.\nCommandline argument missing.\n"\
              "Please enter a positive integer for the number of equations"\
              " to calculate\n\tand run the program again.\n" << std::endl;
  } else if (number > max_number_size) { // arbitrary max number size
    // Too large a number can cause a stack overflow when
    //  equation struct array is created.
    number = max_number_size;
    std::cout << "Too many equation calculations requested.\n"\
              "Only " << max_number_size << " calculations will be "\
              "performed" << std::endl;
  }
  
  return number; // postive int or 0
}

// Assigns address of new arrays to coeffs and roots pointers
void equation_init(equation_t& eq){
  // These arrays are allocated from the heap and exist until deleted.
  // Next version I plan to handle that with a class destructor.
  eq.coeffs = new double[eq.coeffsCount];  
  eq.roots = new double[eq.rootsCount];
  return;
}

// User input is requested and validated
void coeffInput(double& d, char c){
  std::cout << "\nEnter coefficient " << c << ": ";
  std::cin >> d;
  
  // the program wasn't able to assign the user input as a double
  while (!std::cin){
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    std::cout << "\nError! Please enter a valid number for the coefficient.\n";
    std::cout << "\nEnter coefficient " << c << ": ";
    std::cin >> d;
    }

  // clear the input buffer of everything after the double is assigned
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
}

// Operator inputs coefficients.  If zero is entered for coeffiecient a, an
// error message is displayed requesting a new entry.
void readCoeffs(equation_t& eq){

  while (true){  // Runs ad-infinitum until break condition is met.
    coeffInput(*eq.coeffs, 'a');
    if (*eq.coeffs) {break;}  // a must not equal zero 
    else {  // operator entered zero for coefficient a
      std::cout << "\nInvalid entry. Please enter a non-zero "\
           "value for a." << std::endl;
      coeffInput(*eq.coeffs, 'a');
    }
  }

  coeffInput(*(eq.coeffs + 1), 'b');
  coeffInput(*(eq.coeffs + 2), 'c');
  return; 
}

// Computes and returns the discriminant.
double discr(double* coeffs){
  return (coeffs[1]*coeffs[1]-4*coeffs[0]*coeffs[2]);
}

// Gets the discriminant and if it's greater than or equal to 0 
// computes the roots and returns true.
void equSolver(equation_t& eq){
  double discriminant = discr(eq.coeffs);
  
  if (discriminant >= 0){
    eq.roots[0] = (-eq.coeffs[1] + sqrt(discriminant))/(2*eq.coeffs[0]);
    eq.roots[1] = (-eq.coeffs[1] - sqrt(discriminant))/(2*eq.coeffs[0]);
  }
 
  eq.rootsExist = (discriminant >= 0)?true:false;
  return; 
}

void outResults(equation_t& eq, std::ofstream& outStream){
  if (eq.rootsExist){
    // Results are appended to the opened file.
    outStream << "\nQuadratic equation with the following coefficients:"
              << std::endl << "a: " << eq.coeffs[0] << "; b: " << eq.coeffs[1] 
              << "; c: " << eq.coeffs[2] << std::endl
              << "has the following roots" << std::endl << "Root1: " 
              << eq.roots[0] << "; Root2: " << eq.roots[1] << ";\n"
              << std::endl;
  } else {
    std::cout << "\nQuadratic equation with the following coefficients:" 
              << std::endl << "a: " <<eq.coeffs[0] << "; b: " << eq.coeffs[1]
              << "; c: " << eq.coeffs[2]  << std::endl 
              << "has no roots in the real domain.\n" << std::endl;
  }
  return;
}

// Returns resources to the heap in case the OS doesn't automatically 
//  retreive them.
void resource_cleanup(equation_t* eq, int number){
  for (int i=0; i<number;i++){
    
    delete[] (eq+i)->coeffs; 
    delete[] (eq+i)->roots;
    
    // not really needed but good practice i think
    (eq+i)->coeffs = nullptr;
    (eq+i)->roots = nullptr;
  }
}
