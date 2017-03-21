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

/*TODO 1:  validate inputs. Can either use template gist or notebook loop tests.
  :: http://www.cplusplus.com/forum/beginner/13044/#msg62827
   TODO 2:  revise comments based on my google doc's guidelines
 */
 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <stdlib.h>

// https://isocpp.org/wiki/faq/templates
//----------------------------------------------------------------------------
// This is the type definition for an input_t type
template <typename Tom> 
struct input_t
  {
  Tom& n; // input_t can modify by reference
  // constructor: implicit type conversion is not allowed
  explicit input_t( Tom& n ): n( n ) { } 
  input_t( const input_t <Tom> & i ): n( i.n ) { }
  };

//----------------------------------------------------------------------------
// inline is a compiler directive that looks for every instance of input(T& x) and 
//  replace it with the code: input_t <T> result(x); return result;
//  **I think this might be equivalent to a lambda expression**
template <typename Tim>
inline
input_t <Tim>
input( Tim& n )
  {
  input_t <Tim> result( n );
  return result;
  }

//----------------------------------------------------------------------------
// we're overriding the istream >> operator to perform input validation
// lvalue is the istream object (cin), rvalue is an input_t type object 
//   initialized with the user input.
// i is a struct of type input_t which is cast with type T
template <typename Tad> 
std::istream& operator >> ( std::istream& ins, const input_t <Tad> & i )
  {
  // Read a line (terminated by ENTER|NEWLINE) from the user
  std::string s;
  std::getline( ins, s );

  // Get rid of any trailing whitespace
  s.erase( s.find_last_not_of( " \f\n\r\t\v" ) + 1 );

  // Read it into the target type
  std::istringstream ss( s ); // constructor passed the string variable
  ss >> i.n;

  // Check to see that there is nothing left over
  if (!ss.eof())
    // something was left over.
    ins.setstate( std::ios::failbit );

  return ins;
  }

//----------------------------------------------------------------------------
// defaults to a quadratic 
struct equation_t {
  int coeffsCount = 3;  // Override Counts for other degree polynomials
  int rootsCount = 2;  
  double* coeffs;  // pointer to coefficient array
  double* roots;  // pointer to roots array
  bool rootsExist = false;  // if true then real roots exist
};

// requests user to enter doubles
void coeffInput(double&, char);

// Creates new[] equation arrays for coeffs and roots.
void equation_init(equation_t&);

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

// equation_cleanup releases memory resourses before program exits
void resource_cleanup(equation_t*, int);


int main(int argc, char* argv[]) {
  // local constants
  const char* OUTPUT_FILE = "results.dat";

  // number of equations to calculate
  int number = std::abs((argc > 1)?std::atoi(argv[1]):0);  
  
  // struct array to hold all the equations calculated
  equation_t quadratic[number];
  
  std::ofstream outStream;
  outStream.open(OUTPUT_FILE);  // if file exists, it is overwritten
  
  // Each iteration calculates one quadratic equation.
  for (int i=0; i < number; i++){
    // initialize this quadratic equation
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
  if (number > 1){
    std::cout << "( " << number << " ) equation" << ((number == 1)?" ":"s ")
                  << "calculated. Have a nice day!" << std::endl;
  } else {
    std::cout << "Invalid commandline argument. Please enter a positive "\
               "integer for the number of equations to calculate." << std::endl;
  }
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

void coeffInput(double& d, char c){
  std::cout << "\nEnter coefficient " << c << ": ";
  std::cin >> d;
  while (!std::cin){
    std::cin.clear();
    std::cout << "Error! Please enter a valid number for the coefficient.\n";
    std::cout << "\nEnter coefficient " << c << ": ";
    std::cin >> d;
  }
  
}

// Operator inputs coefficients.  If zero is entered for coeffiecient a, an
// error message is displayed requesting a new entry.
void readCoeffs(equation_t& eq){

  while (true){  // Runs ad-infinitum until break condition is met.
     
    coeffInput(*eq.coeffs, 'a');
    if (*eq.coeffs) break;  // a must not equal zero 
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
    outStream << "Quadratic equation with the following coefficients:"
              << std::endl << "a: " << eq.coeffs[0] << "; b: " << eq.coeffs[1] 
              << "; c: " << eq.coeffs[2] << std::endl
              << "has the following roots" << std::endl << "Root1: " 
              << eq.roots[0] << "; Root2: " << eq.roots[1] << ";" << std::endl 
              << std::endl;
  } else {
    std::cout << std::endl << "Quadratic equation with the following coefficients:" 
              << std::endl << "a: " <<eq.coeffs[0] << "; b: " << eq.coeffs[1]  
              << "; c: " << eq.coeffs[2]  << std::endl  
              << "has no roots in the real domain." << std::endl << std::endl;
  }
  return;
}

// comments needed?
void resource_cleanup(equation_t* eq, int number){
  for (int i=0; i<number;i++){
    
    delete[] (eq+i)->coeffs; 
    delete[] (eq+i)->roots;
    
    // not really needed but good practice i think
    (eq+i)->coeffs = nullptr;
    (eq+i)->roots = nullptr;
  }
}
