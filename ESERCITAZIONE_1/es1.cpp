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

double error(int N, double sum, double sum2)
{
   if (N == 0)
   {
      return 0;
   }
   else
   {
      return sqrt(1. / N * ((1. / N) * sum2 - pow(((1. / N) * sum), 2)));
   }
}

double chi2(int n_i, int n, int M)
{
   double y = n / M;
   return pow(n_i - y, 2) / y;
}

int place_in_interval(double r, int M)
{

   int result = 0;

   for (int i = 0; i < M; i++)
   {
      if (r >= i * 1.0 / M && r < (i + 1) * 1.0 / M)
      {
         result = i;
         break;
      }
   }
   return result;
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

   int M = 100000; // lanci totali
   int N = 100;    // numero blocchi
   int L = M / N;  // L lanci per blocco
   ofstream out("Filephy.txt");
   ofstream vout("Filephy_var.txt");
   ofstream wout("chi2_out.txt");

   double A = 0;
   double A2 = 0;
   double mean = 0;
   double mean2 = 0;
   double var = 0;
   double var2 = 0;
   double err_mean;
   double err_var;

   for (int i = 1; i <= N; i++)
   {
      double ran;
      err_mean = 0;
      err_var = 0;
      A = 0;
      A2 = 0;
      for (int j = 1; j <= M / N; j++)
      {
         ran = rnd.Rannyu();
         A = A + ran;
         A2 = A2 + ran * ran;
      }
      mean = mean + A / L;
      mean2 = mean2 + A / L * A / L;
      err_mean = error(i, mean, mean2);

      var = var + A2 / L - pow(A / L, 2);
      var2 = var2 + pow(A2 / L - pow(A / L, 2), 2);
      err_var = error(i, var, var2);

      out << i << "\t" << mean / i - 0.5 << "\t" << err_mean << endl;
      vout << i << "\t" << (var) / i - 1. / 12. << "\t" << err_var << endl;
   }

   // chi^2part
   int m = 100;
   int n = pow(10, 4);
   vector<int> X(m); // X[i]==n_i
   double chi2_sum;

   for (int j = 1; j <= 100; j++)
   {
      chi2_sum = 0;
      fill(X.begin(), X.end(), 0);

      for (int i = 1; i <= n; i++)
      {
         double ran;
         ran = rnd.Rannyu();
         int p = place_in_interval(ran, m);
         X[p]++;
      }

      for (int i = 0; i < m; i++)
      {
         chi2_sum = chi2_sum + chi2(X[i], n, m);
      }

      wout << j << "\t" << chi2_sum << endl;
   }

   out.close();
   vout.close();
   wout.close();

   // SVUTO E COMPRIMO IL VECTOR
   X.clear();
   X.shrink_to_fit();

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
