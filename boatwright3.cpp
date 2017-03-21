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
#include <cmath>
#include <fstream>
#include <stdlib.h>
//----------------------------------------------------------------------------
template <typename T>
struct input_t
  {
  mutable T& n;
  explicit input_t( T& n ): n( n ) { }
  input_t( const input_t <T> & i ): n( i.n ) { }
  };

//----------------------------------------------------------------------------
template <typename T>
inline
input_t <T>
input( T& n )
  {
  input_t <T> result( n );
  return result;
  }

//----------------------------------------------------------------------------
template <typename T>
istream& operator >> ( istream& ins, const input_t <T> & i )
  {
  // Read a line (terminated by ENTER|NEWLINE) from the user
  string s;
  getline( ins, s );

  // Get rid of any trailing whitespace
  s.erase( s.find_last_not_of( " \f\n\r\t\v" ) + 1 );

  // Read it into the target type
  istringstream ss( s );
  ss >> i.n;

  // Check to see that there is nothing left over
  if (!ss.eof())
    ins.setstate( ios::failbit );

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

  // local variables
  int number = (argc > 1)?std::atoi(argv[1]):0;  // number of equations to calculate
  
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
    
    // Directs output to either the file or stdout respectively, based on rootsExist.
    outResults(quadratic[i], outStream);
  }
  outStream.close();
  resource_cleanup(quadratic, number);
  std::cout << "( " << number << " ) equation" << ((number == 1)?" ":"s ")
                  << "calculated. Have a nice day!" << std::endl;
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
void readCoeffs(equation_t& eq){

  while (true){  // Runs ad-infinitum until break condition is met.
    std::cout << "\nEnter coefficient a: "; 
    std::cin >> *eq.coeffs;
    if (*eq.coeffs) break;  // a must not equal zero 
    else {  // operator entered zero for coefficient a
      std::cout << "\nInvalid entry. Please enter a non-zero "\
           "value for a." << std::endl;
    }
  }
  std::cout << "\nEnter coefficient b: "; 
  std::cin >> *(eq.coeffs + 1);
  std::cout << "\nEnter coefficient c: "; 
  std::cin >> *(eq.coeffs + 2);
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
                        << "has the following roots" << std::endl << "Root1: " << eq.roots[0] 
                        << "; Root2: " << eq.roots[1] << ";" << std::endl << std::endl;
  } else {
    std::cout << std::endl << "Quadratic equation with the following coefficients:" 
                    << std::endl << "a: " <<eq.coeffs[0] << "; b: " << eq.coeffs[1] << "; c: " 
                   << eq.coeffs[2]  << std::endl << "has no roots in the real domain." 
                    << std::endl << std::endl;
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
