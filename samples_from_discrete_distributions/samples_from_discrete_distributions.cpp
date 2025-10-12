// A simple C++ program to generate random samples from the geometric distribution and the zeta distribution
// Version 0.1. December 2013
// Author: Ramon Ferrer-i-Cancho

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits> 

using namespace std;

const int length = 1000;

inline double random_uniform_deviate_between_0_and_1() {
   return random()/double(RAND_MAX);
}

// Inspired by Dagpunar's book (?)
class geometric_deviate {
private:
   double b;
public: 
   geometric_deviate(double pi) {
      assert(0 < pi);
      assert(pi < 1);
      b = log(1 - pi);
   }
   long int generate() const {
      double U,X;
      do {
         do { 
            U = random_uniform_deviate_between_0_and_1();
         } while (U == 0);
         X = floor(log(U)/b + 1);
         assert(X >= 1);
      } 
      while (X > numeric_limits<long int>::max());
      return (long int)X;
   }
};

// Based on Devroye's famous book, Section X.6. "The Zipf distribution". pp. 550.
class zeta_deviate {
private:
   double a_minus_1;
   double b;
   double c;
public: 
   zeta_deviate(double a) {
      assert(a > 1);
      a_minus_1 = a - 1;
      b = pow(2, a - 1);
      c = -1/(a_minus_1);
   }
   long int generate() const {
      double U,V,X,T;
      do {
         U = random_uniform_deviate_between_0_and_1();
         V = random_uniform_deviate_between_0_and_1();
         X = floor(pow(U, c));
         T = pow(1 + 1/X, a_minus_1);
      } 
      while (V*X*(T-1)/(b - 1) > T/b or X > numeric_limits<long int>::max());
      return (long int)X;
   }
};

template <class T> void generate_sample(const string &name, const T &d, double parameter) {
   ostringstream outs;
   outs << "./data/sample_of_" << name << "_with_parameter_" << parameter << ".txt";
   ofstream out(outs.str().c_str());
   assert(out);
   cout << "Producing " << outs.str() << endl;
   int sum = 0;
   for (int i = 0; i < length; ++i) {
       long int x = d.generate();
       assert(x >= 1);
       out << x << endl; 
       sum += x;
   }  
   cout << "   Mean value of samples: " << double(sum)/length << " " << 1/parameter << endl;
}
 
void generate_sample_zeta(double parameter);
void generate_sample_geometric(double parameter);

int main() {
   generate_sample_geometric(0.05);
   generate_sample_geometric(0.1);
   generate_sample_geometric(0.2);
   generate_sample_geometric(0.4);
   generate_sample_geometric(0.8);
   generate_sample_zeta(1.5);
   generate_sample_zeta(2);
   generate_sample_zeta(2.5);
   generate_sample_zeta(3);
   generate_sample_zeta(3.5);
}

void generate_sample_geometric(double parameter) {
   geometric_deviate d(parameter);
   generate_sample("geometric", d, parameter);
}

void generate_sample_zeta(double parameter) {
   zeta_deviate d(parameter);
   generate_sample("zeta", d, parameter);
}

