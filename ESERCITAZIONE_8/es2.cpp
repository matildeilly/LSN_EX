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
{

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

double equilibration(double mu, double sigma, double &x_old, Random *rnd, double &larghezza)
{
   int nsteps = 1000;
   int test_each = 50;
   int acceptance = 0;
   int count = 0;
   double epsilon = 0.05;

   for (int i = 0; i < nsteps; i++)
   {
      acceptance += test_acceptance(mu, sigma, x_old, rnd, larghezza);
      count++;
      if (count == test_each)
      {
         if (double(acceptance) / count > 0.6)
         {
            larghezza = larghezza * (1 + epsilon); // se l'accettazione è grande è perchè mi muovo poco --> aumento la larghezza
         }
         if (double(acceptance) / count < 0.4)
         {
            larghezza = larghezza * (1 - epsilon); // per aumentare l'accettazione mi devo muovere di meno --> diminuisco la larghezza
         }
         // cout<<double(acceptance)/count<<endl;
         count = 0;
         acceptance = 0;
      }
   }
   return larghezza;
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

double energy(int nsteps, double &x_in, double mu, double sigma, Random *rnd, double &length_x)
{ // es1 in forma di funzione

   length_x = equilibration(mu, sigma, x_in, rnd, length_x); // EQUILIBRAZIONE (x_length and x_in change)
   double x;
   double energy_sum = 0;

   for (int i = 0; i < nsteps; i++)
   {

      if (i == 0)
      {
         x = x_in;
      }
      if (i != 0)
      {
         x = avanza_linear(mu, sigma, x, rnd, length_x);
      }

      energy_sum += V(x) + kin(mu, sigma, x);
      x_in = x;
   }
   return energy_sum / double(nsteps);
}

void final(int nsteps, double &x_in, double mu, double sigma, Random *rnd, double length_x, ofstream &x_out, ofstream &eminout)
{ // una volta trovati i mu e sigma che minimizzano l'energia, calcolo l'energia media e la varianza (con previa equilibrazione) 

   length_x = equilibration(mu, sigma, x_in, rnd, length_x); // EQUILIBRAZIONE (x_length and x_in change)
   double x;
   double energy_sum = 0;
   int N = 100;
   int L = nsteps / N;
   int pos_in_vec = 0;
   vector<double> E(N);

   for (int i = 0; i < nsteps; i++)
   {

      if (i == 0)
      {
         x = x_in;
      }
      if (i != 0)
      {
         x = avanza_linear(mu, sigma, x, rnd, length_x);
      }
      x_out << x << endl;
      energy_sum += V(x) + kin(mu, sigma, x);
      x_in = x;

      if ((i + 1) % L == 0)
      {
         E[pos_in_vec] = energy_sum / L;
         eminout << pos_in_vec + 1 << "\t" << Vector_mean(E, pos_in_vec + 1) << "\t" << Vector_variance(E, pos_in_vec + 1) << endl;
         energy_sum = 0.;
         pos_in_vec++;
      }
   }
   E.clear();
   E.shrink_to_fit();
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

   ofstream eout("SA_energy.txt");
   ofstream msout("mu_sigma.txt");
   ofstream xout("x_distribution.txt");
   ofstream eminout("energy_min.txt");

   double mu = 1.1;
   double sigma = 0.4;
   double mu_sigma_length = 0.05;   //larghezza di campionamento di mu e sigma
   double x_in = 1.1;
   double x_length = 0.8;

   double T_0 = 10;
   double T_min = 0.0001;
   double rate = 0.99;
   int steps_per_temp = 100000;
   double mu_proposed, sigma_proposed;
   int index = 0;
   double E_old = energy(steps_per_temp, x_in, mu, sigma, &rnd, x_length);
   double E_new;
   double beta;

   for (double T = T_0; T >= T_min; T = T * rate)
   {
      beta = 1.0 / T;
      mu_proposed = rnd.Rannyu(mu - mu_sigma_length / 2., mu + mu_sigma_length / 2.);
      sigma_proposed = rnd.Rannyu(sigma - mu_sigma_length / 2., sigma + mu_sigma_length / 2.);
      E_new = energy(steps_per_temp, x_in, mu_proposed, sigma_proposed, &rnd, x_length);
      double A = min(1., exp(-beta * (E_new - E_old)));
      double soglia = rnd.Rannyu();
      if (soglia <= A)
      {
         E_old = E_new;
         mu = mu_proposed;
         sigma = sigma_proposed;
      }
      eout << index << "\t" << E_new << endl;
      msout << mu << "\t" << sigma << endl;
      index++;
   }

   final(steps_per_temp, x_in, mu, sigma, &rnd, x_length, xout, eminout);

   eout.close();
   msout.close();
   xout.close();
   eminout.close();

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