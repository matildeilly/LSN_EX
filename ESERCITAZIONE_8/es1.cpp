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

double psi(double mu, double sigma, double x)
{
   return exp(-(x - mu) * (x - mu) / (2 * sigma * sigma)) + exp(-(x + mu) * (x + mu) / (2 * sigma * sigma));
}

double V(double x)
{
   return pow(x, 4) - 2.5 * x * x;
}

double kin(double mu, double sigma, double x)
{
   double num = pow((x - mu) / (sigma * sigma), 2) * exp(-pow((x - mu), 2) / (2 * sigma * sigma)) + pow((x + mu) / (sigma * sigma), 2) * exp(-pow((x + mu), 2) / (2 * sigma * sigma));
   return -0.5 * (num / psi(mu, sigma, x) - pow(sigma, -2));
}

double min(double x, double y)
{
   if (x < y)
      return x;
   else
      return y;
}

double avanza_linear(double mu, double sigma, double x_old, Random *rnd, double larghezza)
{

   double x_new;
   x_new = rnd->Rannyu(x_old - larghezza / 2., x_old + larghezza / 2.);

   double A = min(1., pow(psi(mu, sigma, x_new), 2) / pow(psi(mu, sigma, x_old), 2));

   double soglia = rnd->Rannyu();
   if (soglia <= A)
   {
      return x_new;
   }

   else
   {
      return x_old;
   }
}

int test_acceptance(double mu, double sigma, double &x_old, Random *rnd, double larghezza)
{ // serve per verificare che l'accettazione sia circa del 50 %

   double x_new;
   x_new = rnd->Rannyu(x_old - larghezza / 2., x_old + larghezza / 2.);

   double A = min(1., pow(psi(mu, sigma, x_new), 2) / pow(psi(mu, sigma, x_old), 2));

   double soglia = rnd->Rannyu();
   if (soglia <= A)
   {
      x_old = x_new;
      return 1;
   }

   else
   {
      return 0;
   }
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

   ofstream out("metro_energy.txt");
   ofstream xout("metro_x.txt");

   int M = pow(10, 5);
   int N = 100;
   int L = M / N;

   double mu = 1.1;
   double sigma = 0.4;

   double length = 1.6; // parametro del metropolis

   double x;
   double x_in = 1.1;
   double energy_sum = 0;

   vector<double> energy(N);
   int pos_in_vec = 0;
   // int test = 0;

   for (int i = 0; i < M; i++)
   {

      if (i == 0)
      {
         x = x_in;
      }

      if (i != 0)
      {
         x = avanza_linear(mu, sigma, x, &rnd, length);
         // test += test_acceptance(mu,sigma,x,&rnd,length);
         // cout<<x<<endl;
      }

      xout << x << endl;
      energy_sum += V(x) + kin(mu, sigma, x);

      if ((i + 1) % L == 0)
      {
         energy[pos_in_vec] = energy_sum / L;
         out << pos_in_vec + 1 << "\t" << Vector_mean(energy, pos_in_vec + 1) << "\t" << Vector_variance(energy, pos_in_vec + 1) << endl;
         energy_sum = 0.;
         pos_in_vec++;

         // cout<<double(test)/L<<endl; // * * per verificare l'accettazione
         // test = 0;
      }
   }

   out.close();
   xout.close();

   // SVUTO E COMPRIMO I VECTORS
   energy.clear();
   energy.shrink_to_fit();

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
