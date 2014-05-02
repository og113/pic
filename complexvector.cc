// Illustration of compiler <complex> class.
#include <iostream>
#include <complex>
#include <vector>
using namespace std;

int main(){

  complex<double> a(-5.,1.),b,i(0.0,1.0);
  vector <complex<double>> v(2);

  cout << "Enter b: ";
  cin >> b;

  cout << "a = " << a << "\n";
  cout << "b = " << b << "\n";

  cout << "a + b = " << a + b << "\n";
  cout << "a * b = " << a * b << "\n";
  cout << "a / b = " << a / b << "\n";
  cout << "|a| = "   << abs(a) << "\n";
  cout << "complex conjugate of a = " << conj(a) << "\n";
  cout << "norm of a = " << norm(a) << "\n";
  cout << "abs of a = " << abs(a) << "\n";
  cout << "exp(a) = " << exp(a) << "\n";
  
  v[0] = a;
  v[1] = b;
  cout << "v = " << "[" << v[0] << ";" << v[1] << "]" << endl;
  
}
