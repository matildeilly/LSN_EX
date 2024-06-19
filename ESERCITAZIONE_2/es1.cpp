/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <vector>

using namespace std;

double function_integral(double x)
{
   return (M_PI / 2 * cos(M_PI * x / 2));
}

double p_function(double y)
{ // p(x) indica la distribuzione di probabilità da campionare; questa funziona rappresenta l'inversa della cumulativa di p(x)
   return 1 - sqrt(1 - y);
}

double g_function(double x)
{ // f(x)/p(x); la valuto in x estratta secondo la distribuzione di probabilità p(x)
   return M_PI / 4 * cos(M_PI / 2 * x) / (1 - x);
}

double Vector_mean(vector<double> V, int dim)
{
   double sum = 0;
   for (int i = 0; i < dim; i++)
   {
      sum = sum + V[i];
   }
   return sum / dim;
}

double Vector_variance(vector<double> V, int dim)
{
   if (dim == 1)
   {
      return 0;
   }
   else
   {
      double sum = 0;
      double sum2 = 0;
      for (int i = 0; i < dim; i++)
      {
         sum = sum + V[i];
         sum2 = sum2 + pow(V[i], 2);
      }

      return sqrt(1.0 / dim * (1.0 / dim * sum2 - pow(1.0 / dim * sum, 2)));
   }
}

int main(int argc, char *argv[])
{

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open())
   {
      Primes >> p1 >> p2;
   }
   else
      cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open())
   {
      while (!input.eof())
      {
         input >> property;
         if (property == "RANDOMSEED")
         {
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed, p1, p2);
         }
      }
      input.close();
   }
   else
      cerr << "PROBLEM: Unable to open seed.in" << endl;

   rnd.SaveSeed();

   int M = pow(10, 4); // numero di throws
   int N = 100;        // numero di blocchi
   int L = M / N;      // numero di throws per blocco

   vector<double> I(N);
   vector<double> I2(N);

   ofstream out("Integral_unif_out.txt");
   ofstream out2("Integral_importance_out.txt");

   for (int i = 0; i < N; i++)
   {
      double ran;
      double ran2;
      double sum = 0;
      double sum2 = 0;

      for (int j = 1; j <= L; j++)
      {
         ran = rnd.Rannyu();
         ran2 = p_function(rnd.Rannyu());
         sum += function_integral(ran);
         sum2 += g_function(ran2);
      }
      I[i] = sum / L;
      I2[i] = sum2 / L;
      out << i + 1 << "\t" << Vector_mean(I, i + 1) << "\t" << Vector_variance(I, i + 1) << endl;
      out2 << i + 1 << "\t" << Vector_mean(I2, i + 1) << "\t" << Vector_variance(I2, i + 1) << endl;
   }

   out.close();
   out2.close();

   // SVUTO E COMPRIMO I VECTORS
   I.clear();
   I.shrink_to_fit();
   I2.clear();
   I2.shrink_to_fit();
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
