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

struct position
{
   double x, y, z;

   position() : x(0), y(0), z(0) {}

   position(double x1, double x2, double x3)
   {
      x = x1;
      y = x2;
      z = x3;
   }

   double dist()
   {
      return sqrt(x * x + y * y + z * z);
   }
};

position avanza(position A, double ran)
{

   if (ran >= 1 && ran < 2)
   {
      return position(A.x - 1, A.y, A.z);
   }
   if (ran >= 2 && ran < 3)
   {
      return position(A.x + 1, A.y, A.z);
   }
   if (ran >= 3 && ran < 4)
   {
      return position(A.x, A.y - 1, A.z);
   }
   if (ran >= 4 && ran < 5)
   {
      return position(A.x, A.y + 1, A.z);
   }
   if (ran >= 5 && ran < 6)
   {
      return position(A.x, A.y, A.z - 1);
   }
   else
   {
      return position(A.x, A.y, A.z + 1);
   }
}

position avanza_continous(position A, Random *rnd)
{

   double theta = rnd->Rannyu(0., 2 * M_PI);
   double phi = acos(2 * rnd->Rannyu() - 1.);

   double x = sin(phi) * cos(theta);
   double y = sin(phi) * sin(theta);
   double z = cos(phi);

   return position(A.x + x, A.y + y, A.z + z);
}

double Vector_mean(vector<double> V)
{
   double sum = 0;
   for (int i = 0; i < V.size(); i++)
   {
      sum = sum + V[i];
   }
   return sum / V.size();
}

double Vector_variance(vector<double> V)
{
   double sum = 0;
   double sum2 = 0;
   double dim = V.size();
   for (int i = 0; i < dim; i++)
   {
      sum += V[i];
      sum2 += pow(V[i], 2);
   }

   return sqrt(1.0 / dim * (1.0 / dim * sum2 - pow(1.0 / dim * sum, 2)));
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

   int M = pow(10, 4); // numero di volte in cui ripeto il random walk
   int N_blocks = 100;
   int L = M / N_blocks;

   ofstream out("Discreto_out.txt");
   ofstream Cout("Continuo_out.txt");

   vector<position> V_i(M);
   vector<position> V_previous(M);
   vector<position> V_i_cont(M);
   vector<position> V_previous_cont(M);

   for (int i = 0; i <= 100; i++) // numero di passi totali
   {

      if (i == 0)
      {
         for (int j = 0; j < M; j++)
         {
            V_i[j].x = 0;
            V_i[j].y = 0;
            V_i[j].z = 0;
            V_i_cont[j].x = 0;
            V_i_cont[j].y = 0;
            V_i_cont[j].z = 0;
         }
         out << i << "\t" << 0 << "\t" << 0 << endl;
         Cout << i << "\t" << 0 << "\t" << 0 << endl;
      }
      else
      {

         double sum_dist = 0;
         vector<double> D(N_blocks);
         int pos_in_vec = 0;
         double sum_dist_cont = 0;
         vector<double> D_cont(N_blocks);
         int pos_in_vec_cont = 0;

         for (int j = 1; j <= M; j++)
         {

            V_i[j - 1] = avanza(V_previous[j - 1], rnd.Rannyu(1, 7));
            V_i_cont[j - 1] = avanza_continous(V_previous_cont[j - 1], &rnd);
            sum_dist += pow(V_i[j - 1].dist(), 2);
            sum_dist_cont += pow(V_i_cont[j - 1].dist(), 2);

            if (j % L == 0)
            {
               D[pos_in_vec] = sqrt(sum_dist / L);
               pos_in_vec++;
               sum_dist = 0;
               D_cont[pos_in_vec_cont] = sqrt(sum_dist_cont / L);
               pos_in_vec_cont++;
               sum_dist_cont = 0;
            }
         }
         out << i << "\t" << Vector_mean(D) << "\t" << Vector_variance(D) << endl;
         Cout << i << "\t" << Vector_mean(D_cont) << "\t" << Vector_variance(D_cont) << endl;
      }

      V_previous = V_i;
      V_previous_cont = V_i_cont;
   }

   out.close();
   Cout.close();

   //   SVUTO E COMPRIMO I VECTORS
   V_i.clear();
   V_i.shrink_to_fit();
   V_i_cont.clear();
   V_i_cont.shrink_to_fit();
   V_previous.clear();
   V_previous.shrink_to_fit();
   V_previous_cont.clear();
   V_previous_cont.shrink_to_fit();

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
