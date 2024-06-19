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

double max(double x, double y)
{
   if (x > y)
   {
      return x;
   }
   else
   {
      return y;
   }
}

double C_dir_value(double K, double r, double T, double S_0, double sigma, double W /*W gaussiana N(0,T)*/)
{

   double S_T = S_0 * exp((r - sigma * sigma / 2.) * T + sigma * W);

   return max(0, S_T - K) * exp(-r * T);
}

double P_dir_value(double K, double r, double T, double S_0, double sigma, double W /*W gaussiana N(0,T)*/)
{

   double S_T = S_0 * exp((r - sigma * sigma / 2.) * T + sigma * W);

   return max(0, K - S_T) * exp(-r * T);
}

double C_interv_value(double K, double r, double T, double S_0, double sigma, Random *rnd)
{

   double delta_t = T / 100.; // poichè l'intervallo [0,T] è diviso in 100 intervalli uguali
   double A = 1.;
   double Z;

   for (int i = 1; i <= 100; i++)
   {
      Z = rnd->Gauss(0., 1.);
      A = A * exp((r - sigma * sigma / 2.) * delta_t + sigma * Z * sqrt(delta_t));
   }

   double S_T = S_0 * A;

   return max(0, S_T - K) * exp(-r * T);
}

double P_interv_value(double K, double r, double T, double S_0, double sigma, Random *rnd)
{

   double delta_t = T / 100.; // poichè l'intervallo [0,T] è diviso in 100 intervalli uguali
   double A = 1.;
   double Z;

   for (int i = 1; i <= 100; i++)
   {
      Z = rnd->Gauss(0., 1.);
      A = A * exp((r - sigma * sigma / 2.) * delta_t + sigma * Z * sqrt(delta_t));
   }

   double S_T = S_0 * A;

   return max(0, K - S_T) * exp(-r * T);
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

   int M = pow(10, 6); // numero di tiri
   int N = 100;        // numero blocchi
   int L = M / N;      // tiri per blocco

   double S_0 = 100.;
   double T = 1.;
   double K = 100.;
   double r = 0.1;
   double sigma = 0.25;
   double W;

   vector<double> C_dir(N);
   vector<double> P_dir(N);
   vector<double> C_interv(N);
   vector<double> P_interv(N);

   ofstream C1out("C_dir_out.txt");
   ofstream P1out("P_dir_out.txt");
   ofstream C2out("C_interv_out.txt");
   ofstream P2out("P_interv_out.txt");

   for (int i = 0; i < N; i++)
   {

      double sum_C_dir = 0.;
      double sum_P_dir = 0.;
      double sum_C_interv = 0.;
      double sum_P_interv = 0.;

      for (int j = 1; j <= L; j++)
      {

         W = rnd.Gauss(0, T);

         sum_C_dir += C_dir_value(K, r, T, S_0, sigma, W);
         sum_P_dir += P_dir_value(K, r, T, S_0, sigma, W);
         sum_C_interv += C_interv_value(K, r, T, S_0, sigma, &rnd);
         sum_P_interv += P_interv_value(K, r, T, S_0, sigma, &rnd);
      }
      C_dir[i] = sum_C_dir / L;
      P_dir[i] = sum_P_dir / L;
      C_interv[i] = sum_C_interv / L;
      P_interv[i] = sum_P_interv / L;

      C1out << i + 1 << "\t" << Vector_mean(C_dir, i + 1) << "\t" << Vector_variance(C_dir, i + 1) << endl;
      P1out << i + 1 << "\t" << Vector_mean(P_dir, i + 1) << "\t" << Vector_variance(P_dir, i + 1) << endl;
      C2out << i + 1 << "\t" << Vector_mean(C_interv, i + 1) << "\t" << Vector_variance(C_interv, i + 1) << endl;
      P2out << i + 1 << "\t" << Vector_mean(P_interv, i + 1) << "\t" << Vector_variance(P_interv, i + 1) << endl;
   }

   C1out.close();
   P1out.close();
   C2out.close();
   P2out.close();

   // SVUTO E COMPRIMO I VECTORS
   C_dir.clear();
   C_dir.shrink_to_fit();
   P_dir.clear();
   P_dir.shrink_to_fit();
   C_interv.clear();
   C_interv.shrink_to_fit();
   P_interv.clear();
   P_interv.shrink_to_fit();

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
