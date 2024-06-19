#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <vector>
#include <climits>

using namespace std;

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

   int M = pow(10, 4);  // numero lanci totali
   int N = 100;         // numero di blocchi
   int L = M / N;       // L lanci per blocco
   vector<double> P(N); // vector contenente la stima di M_PI per ogni blocco

   double x1, x2, x_cos, y_cos, coseno;

   ofstream out("Buffon_out.txt");

   for (int i = 0; i < N; i++)
   {
      int N_hit;
      N_hit = 0;

      for (int j = 1; j <= L; j++)
      {
         x1 = rnd.Rannyu(); // primo estremo della sbarra
         x_cos = rnd.Rannyu(-1, 1);
         y_cos = rnd.Rannyu();
         // ho estratto così un angolo da 0 a pi (non 2pi perchè tanto voglio il coseno)
         coseno = x_cos / sqrt(x_cos * x_cos + y_cos * y_cos);

         x2 = x1 + 0.9 * coseno; // secondo punto della sbarra (la sbarra ha lunghezza 0.9)

         if (x2 >= 1 || x2 <= 0)
         {
            N_hit++;
         } // la griglia ha spaziatura 1
      }

      P[i] = 2 * 0.9 * L / (N_hit * 1);

      out << i + 1 << "\t" << Vector_mean(P, i + 1) << "\t" << Vector_variance(P, i + 1) << endl;
   }

   out.close();

   // SVUTO E COMPRIMO IL VECTOR
   P.clear();
   P.shrink_to_fit();

   return 0;
}
