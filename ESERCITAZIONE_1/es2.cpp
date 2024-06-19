#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>
#include <vector>

using namespace std;

double prob_exponential(double y, double lambda)
{

    return -1 / lambda * log(1 - y);
}

double prob_Lorentz(double y, double gamma)
{

    if (abs(gamma * tan(M_PI * (y - 1.0 / 2))) > 100)
    {
        return 0;
    }

    else
    {
        return gamma * tan(M_PI * (y - 1.0 / 2));
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

    int M = pow(10, 4);
    double sum2 = 0;
    double sum10 = 0;
    double sum100 = 0;
    double sum2_E = 0;
    double sum10_E = 0;
    double sum100_E = 0;
    double sum2_L = 0;
    double sum10_L = 0;
    double sum100_L = 0;
    int pos2 = 0;
    int pos10 = 0;
    int pos100 = 0;

    vector<double> Unif_1(M);
    vector<double> Unif_2(M);
    vector<double> Unif_10(M);
    vector<double> Unif_100(M);
    vector<double> Exp_1(M);
    vector<double> Exp_2(M);
    vector<double> Exp_10(M);
    vector<double> Exp_100(M);
    vector<double> Lorentz_1(M);
    vector<double> Lorentz_2(M);
    vector<double> Lorentz_10(M);
    vector<double> Lorentz_100(M);

    ofstream Uout("Unif_out.txt");
    ofstream Eout("Exp_out.txt");
    ofstream Lout("Lorentz_out.txt");

    for (int i = 1; i <= 100 * M; i++)
    {
        double ran;
        ran = rnd.Rannyu();

        if (i <= M)
        {
            Unif_1[i - 1] = ran;
            Lorentz_1[i - 1] = prob_Lorentz(ran, 1);
            Exp_1[i - 1] = prob_exponential(ran, 1);
        }

        if (i <= 2 * M)     
        {

            if (i % 2 != 0)     //in questo caso mi serve ancora un altro numero per poter fare la media
            {
                sum2 = sum2 + ran;
                sum2_E = sum2_E + prob_exponential(ran, 1);
                sum2_L = sum2_L + prob_Lorentz(ran, 1);
            }
            else                //ora faccio la media
            {
                sum2 = sum2 + ran;
                sum2_E = sum2_E + prob_exponential(ran, 1);
                sum2_L = sum2_L + prob_Lorentz(ran, 1);
                Unif_2[pos2] = sum2 / 2;
                Lorentz_2[pos2] = sum2_L / 2;
                Exp_2[pos2] = sum2_E / 2;
                sum2 = 0;
                sum2_E = 0;
                sum2_L = 0;
                pos2++;
            }
        }

        if (i <= 10 * M)
        {

            if (i % 10 != 0)
            {
                sum10 = sum10 + ran;
                sum10_E = sum10_E + prob_exponential(ran, 1);
                sum10_L = sum10_L + prob_Lorentz(ran, 1);
            }
            else
            {
                sum10 = sum10 + ran;
                sum10_E = sum10_E + prob_exponential(ran, 1);
                sum10_L = sum10_L + prob_Lorentz(ran, 1);
                Unif_10[pos10] = sum10 / 10;
                Lorentz_10[pos10] = sum10_L / 10;
                Exp_10[pos10] = sum10_E / 10;
                sum10 = 0;
                sum10_E = 0;
                sum10_L = 0;
                pos10++;
            }
        }

        if (i <= 100 * M)
        {

            if (i % 100 != 0)
            {
                sum100 = sum100 + ran;
                sum100_E = sum100_E + prob_exponential(ran, 1);
                sum100_L = sum100_L + prob_Lorentz(ran, 1);
            }
            else
            {
                sum100 = sum100 + ran;
                sum100_E = sum100_E + prob_exponential(ran, 1);
                sum100_L = sum100_L + prob_Lorentz(ran, 1);
                Unif_100[pos100] = sum100 / 100;
                Lorentz_100[pos100] = sum100_L / 100;
                Exp_100[pos100] = sum100_E / 100;
                sum100 = 0;
                sum100_E = 0;
                sum100_L = 0;
                pos100++;
            }
        }
    }

    for (int i = 0; i < M; i++)
    {
        Uout << Unif_1[i] << "\t" << Unif_2[i] << "\t" << Unif_10[i] << "\t" << Unif_100[i] << endl;
        Eout << Exp_1[i] << "\t" << Exp_2[i] << "\t" << Exp_10[i] << "\t" << Exp_100[i] << endl;
        Lout << Lorentz_1[i] << "\t" << Lorentz_2[i] << "\t" << Lorentz_10[i] << "\t" << Lorentz_100[i] << endl;
    }

    Uout.close();
    Eout.close();
    Lout.close();

    // SVUTO E COMPRIMO I VECTORS
    Unif_1.clear();
    Unif_1.shrink_to_fit();
    Unif_2.clear();
    Unif_2.shrink_to_fit();
    Unif_10.clear();
    Unif_10.shrink_to_fit();
    Unif_100.clear();
    Unif_100.shrink_to_fit();

    Exp_1.clear();
    Exp_1.shrink_to_fit();
    Exp_2.clear();
    Exp_2.shrink_to_fit();
    Exp_10.clear();
    Exp_10.shrink_to_fit();
    Exp_100.clear();
    Exp_100.shrink_to_fit();

    Lorentz_1.clear();
    Lorentz_1.shrink_to_fit();
    Lorentz_2.clear();
    Lorentz_2.shrink_to_fit();
    Lorentz_10.clear();
    Lorentz_10.shrink_to_fit();
    Lorentz_100.clear();
    Lorentz_100.shrink_to_fit();

    return 0;
}
