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
#include <cassert>
#include <algorithm>
#include <mpi.h>

using namespace std;

vector<int> ordinato(int dim)
{
    vector<int> ordinatoo(dim);
    for (int i = 0; i < dim; i++)
    {
        ordinatoo[i] = i;
    }
    return ordinatoo;
}

bool check(vector<int> V)
{ // da applicare a tutta la generazione: se ne tolgo uno, togline anche un altro selezionato con select inverso (perchè voglio generazione con numero pari (per crossing))
    for (int i = 0; i < V.size(); i++)
    {
        int count = 0;
        for (int j = 0; j < V.size(); j++)
        {
            if (j == 0 && V[j] != 0)
            {
                return false;
            }
            if (V[j] == i)
            {
                count++;
            }
        }
        if (count != 1)
        {
            return false;
        }
    }
    return true;
}

vector<vector<int>> check_gen(vector<vector<int>> &gen)
{
    vector<vector<int>> checked_gen;
    for (int i = 0; i < gen.size(); i++)
    {
        if (check(gen[i]))
            checked_gen.push_back(gen[i]);
        assert(check(gen[i]) && "Invalid path generated after check");
    }
    if (checked_gen.size() % 2 != 0)
    {                                                            // volgio che la dim sia pari (vedi crossing...)
        checked_gen.insert(checked_gen.begin(), checked_gen[0]); // se non è pari, aggiungo il migliore
    }

    return checked_gen;
}

struct position
{
    double x, y;
    position() : x(0), y(0) {}
    position(double x1, double x2)
    {
        x = x1;
        y = x2;
    }
};

double dist(position a, position b)
{
    return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2));
}

double metric(vector<int> A, vector<position> points)
{
    int dim = A.size();
    vector<position> V(dim);
    for (int i = 0; i < dim; i++)
    {
        V[i] = points[A[i]];
    }
    double add = 0.;
    for (int i = 0; i < (dim - 1); i++)
    { //(dim - 1) così all'ultimo step ho la distanza tra il penultimo e l'ultimo elemento
        add += dist(V[i], V[i + 1]);
    }
    add += dist(V[0], V[dim - 1]);
    return add;
}

int identifica_pos_in_vec(double val, vector<double> V)
{
    double dim = V.size();
    for (int i = 0; i < dim; i++)
    {
        if (V[i] == val)
        {
            return i;
        }
    }
    cerr << "ERROR IN IDENTIFICATION OF POSITION IN VECTOR!!" << endl;
    return -1;
}

vector<vector<int>> sorting(vector<vector<int>> &V, vector<position> points)
{
    int dim = V.size(); // V è la generazione: l'insieme di tutti i percorsi
    vector<double> metriche(dim);
    for (int i = 0; i < dim; i++)
    {
        metriche[i] = metric(V[i], points); // associo ad ogni percorso la sua lunghezza
    }
    vector<double> metriche_to_sort = metriche;
    sort(metriche_to_sort.begin(), metriche_to_sort.end());
    vector<vector<int>> new_vec(dim);
    for (int i = 0; i < dim; i++)
    {
        new_vec[i] = V[identifica_pos_in_vec(metriche_to_sort[i], metriche)];
    }
    return new_vec;
}

vector<int> select(vector<vector<int>> V, Random *rnd)
{
    int M = V.size();
    double p = 2;
    int j = (int)(M * pow(rnd->Rannyu(), p));
    return V[j];
}

void mutation_1(vector<int> &V, Random *rnd)
{
    int dim = V.size();
    int index_1 = (int)rnd->Rannyu(1., dim);
    int index_2 = (int)rnd->Rannyu(1., dim);
    swap(V[index_1], V[index_2]);
    assert(check(V) && "Invalid path generated after mutation_1");
}

void mutation_2(vector<int> &V, Random *rnd)
{
    int N = V.size();
    int m = (int)rnd->Rannyu(1, N - 1); // indice da cui partire
    int n = (int)rnd->Rannyu(1, N - m); // numero di città da prendere

    vector<int> temp(n); // temp è un vettore accessorio...

    for (int i = 0; i < n; ++i)
    { // in temp salvo il blocco
        temp[i] = V[m + i];
    }

    for (int i = m + n; i < N; ++i)
    { // sposto gli elementi del blocco successivo al blocco corrente
        V[i - n] = V[i];
    }

    for (int i = 0; i < n; ++i)
    { // riempio la fine del vettore con il blocco originale
        V[N - n + i] = temp[i];
    }
    assert(check(V) && "Invalid path generated after mutation_2");
}

void mutation_3(vector<int> &V, Random *rnd)
{
    int N = V.size();
    int m = (int)rnd->Rannyu(1, N / 2.);         // indice da cui partire
    int n = (int)rnd->Rannyu(1, N / 2. - m + 1); // dimensione del blocco di città da prendere
    int m2 = (int)rnd->Rannyu(m + n, N - n);     // indice da cui far partire il secondo blocco

    for (int i = 0; i < n; i++)
    {
        swap(V[i + m2], V[i + m]);
    }
    assert(check(V) && "Invalid path generated after mutation_3");
}

void mutation_4(vector<int> &V, Random *rnd)
{
    int N = V.size();
    int m = (int)rnd->Rannyu(1, N);         // indice da cui partire
    int n = (int)rnd->Rannyu(1, N - m + 1); // numero di città da prendere

    for (int i = 0; i < n / 2.; i++)
    {
        swap(V[m + i], V[m + n - 1 - i]);
    }
    assert(check(V) && "Invalid path generated after mutation_4");
}

void kaos(vector<vector<int>> &gen, Random *rnd, double prob)
{ // prob è la probabilità di FARE la mutazione (tipo 10%...)
    int dim = gen.size();
    for (int i = 0; i < dim; i++)
    {
        double discriminator = rnd->Rannyu(0, 1);

        if (discriminator < prob)
        {
            if (discriminator < prob / 4.)
            {
                mutation_1(gen[i], rnd);
            }
            else if (discriminator < prob / 2.)
            {
                mutation_2(gen[i], rnd);
            }
            else if (discriminator < prob * 3. / 4.)
            {
                mutation_3(gen[i], rnd);
            }
            else
            {
                mutation_4(gen[i], rnd);
            }
        }
    }
}

vector<vector<int>> initialize_generation(Random *rnd, int dim_gen, int number_of_points, int number_of_permut)
{
    vector<vector<int>> gen(dim_gen);
    vector<int> ord = ordinato(number_of_points);
    for (int i = 0; i < dim_gen; i++)
    {
        gen[i] = ord;
        for (int j = 0; j < number_of_permut; j++)
        {
            mutation_1(gen[i], rnd);
        }
    }
    return gen;
}

bool search(int value, const vector<int> &V)
{ // const per assicurarci che search non modifichi V e & (by reference) per efficienza
    for (int i = 0; i < V.size(); i++)
    {
        if (V[i] == value)
            return true;
    }
    return false;
}

vector<int> update_small(vector<int> &small_parent, const vector<int> &other_parent)
{
    for (int i = 0; i < other_parent.size(); i++)
    {
        if (!search(other_parent[i], small_parent))
        {
            small_parent.push_back(other_parent[i]);
        }
    }
    return small_parent;
}

int identifica_pos_in_gen(vector<int> V, vector<vector<int>> gen)
{
    double dim = gen.size();
    for (int i = 0; i < dim; i++)
    {
        if (gen[i] == V)
        {
            return i;
        }
    }
    cerr << "ERROR IN IDENTIFICATION OF POSITION IN VECTOR!!" << endl;
    return -1;
}

vector<vector<int>> crossing(vector<vector<int>> &gen, Random *rnd, double prob)
{ // prob di FARE il crossing (> 50%)
    int dim_gen = gen.size();
    vector<vector<int>> next_gen;
    int number_of_cross = dim_gen / 2; //!!!!DIM GEN DEVE ESSERE PARI!!!

    for (int i = 0; i < number_of_cross; i++)
    {
        vector<int> parent_1 = select(gen, rnd);
        vector<int> parent_2 = select(gen, rnd);

        while (identifica_pos_in_gen(parent_2, gen) == identifica_pos_in_gen(parent_1, gen))
        { // ci accertiamo che parent_1 e parent_2 siano diversi
            parent_2 = select(gen, rnd);
        }

        double discriminator = rnd->Rannyu(0, 1);

        if (discriminator > prob)
        {
            next_gen.push_back(parent_1);
            next_gen.push_back(parent_2);
        }
        else
        {
            int N = parent_1.size();              // per forza la stessa di parent_2...
            int cut_off = (int)rnd->Rannyu(1, N); // assicurarsi che il taglio non sia all'inizio o alla fine
            vector<int> small_parent_1(parent_1.begin(), parent_1.begin() + cut_off);
            vector<int> small_parent_2(parent_2.begin(), parent_2.begin() + cut_off);

            // completare i figli con le città mancanti nell'ordine in cui appaiono nell'altro genitore
            update_small(small_parent_1, parent_2);
            update_small(small_parent_2, parent_1);

            next_gen.push_back(small_parent_1);
            next_gen.push_back(small_parent_2);
            assert(check(small_parent_1) && "Invalid path generated after crossing");
            assert(check(small_parent_2) && "Invalid path generated after crossing");
        }
    }
    return next_gen;
}

vector<position> grand_tour(string name)
{
    ifstream file1(name);
    int Npoints = 0;
    double a, b;
    while (file1 >> a >> b)
    { // capisco quanti punti ha il file (non riempio subito il vector con push_back perchè è un vector di position!)
        Npoints++;
    }
    file1.close();
    ifstream file(name); // riapro il file e lo rileggo (devo chiuderlo se no non riempie nulla)
    vector<position> points(Npoints);
    int index = 0;
    while (file >> a >> b)
    {
        points[index].x = a;
        points[index].y = b;
        index++;
    }
    file.close();

    return points;
}

void save_path(vector<int> best_path, vector<position> points, ofstream &out)
{
    for (int i = 0; i < best_path.size(); i++)
    {
        int index = best_path[i];
        out << points[index].x << "\t" << points[index].y << endl;
    }
    out << points[0].x << "\t" << points[0].y << endl;
}

void save_length(vector<int> best_path, vector<position> points, ofstream &out, int index)
{
    out << index << "\t" << metric(best_path, points) << endl;
}

int main(int argc, char *argv[])
{
    int size, rank;
    //* INIZIO CALCOLO PARALLELO
    MPI_Init(&argc, &argv);     
    MPI_Comm_size(MPI_COMM_WORLD, &size);   // quanti cores uso
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //chi è ciasuno

    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");

    if (!Primes.is_open())
    {
        cerr << "PROBLEM: Unable to open Primes" << endl;
    }
    Primes >> p1 >> p2;
    for (int i = 0; i < rank + 1; i++)  // per ogni cores devo prendere primes diversi (altrimenti fanno tutti la stessa cosa...)
    {
        Primes >> p1 >> p2;
    }

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
    {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    rnd.SaveSeed();

    ofstream out("italy_path.txt");

    vector<position> points = grand_tour("cap_prov_ita.dat"); // creo il vettore di tutte le città
    int Npoints = points.size();

    int dim_gen = 500;
    if (rank == 0)
        cout << "Sto incominciando l'inizializzazione!" << endl;
    vector<vector<int>> first_generation = initialize_generation(&rnd, dim_gen, Npoints, 400); // creo la prima generazione
    if (rank == 0)
        cout << "Ho finito l'inizializzazione! :)" << endl;
    int Niter = 1000;
    int Nerasmus = 50; // numero di volte in cui i cores si scambiano il miglior percorso
    double prob_cross = 0.7;
    double prob_mut = 0.2;
    vector<vector<int>> generation = sorting(first_generation, points);

    for (int i = 0; i < Niter; i++)
    {
        generation = crossing(generation, &rnd, prob_cross);
        kaos(generation, &rnd, prob_mut);
        generation = sorting(generation, points);

        if ((i + 1) % Nerasmus == 0)
        {
            for (int j = 0; j < size; j++)
            {
                MPI_Status stat1, stat2; // da informazioni sulla comunicazione 
                int tag1 = 1; // identificatore della prima comunicazione
                int tag2 = 2; // identificatore della seconda comunicazione
                int other_rank; // rank del core con cui il core j scambia il miglior percorso
                vector<int> mesg1(Npoints); // vector dove verrà inviato il miglior percorso del core j
                vector<int> mesg2(Npoints); // vector dove verrà ricevuto il miglior percorso del core other_rank

                if (rank == j)
                {
                    mesg1 = generation[0]; // riempio mesg1 con il miglior percorso del core j
                    do
                    {
                        other_rank = (int)rnd.Rannyu(0, size);
                    } while (other_rank == j);  // impongo che other_rank sia diverso da j
                }
                MPI_Bcast(&other_rank /*ciò che j comunica a tutti*/, 1 /*numero di elementi da comunicare*/, MPI_INTEGER /*tipo del messaggio*/, j /*chi comunica*/, MPI_COMM_WORLD /*a chi viene comunicato*/); 
                MPI_Bcast(mesg1.data(), Npoints, MPI_INTEGER, j, MPI_COMM_WORLD);
                if (rank == other_rank)
                {
                    mesg2 = generation[0];
                }

                MPI_Bcast(mesg2.data() /* "data" perchè voglio fare un passaggio per reference e ho un vector*/, Npoints, MPI_INTEGER, other_rank, MPI_COMM_WORLD); // other_rank dice a tutti chi è il suo miglior percorso

                if (rank == j)
                {
                    MPI_Send(mesg1.data(), mesg1.size(), MPI_INTEGER, other_rank, tag1, MPI_COMM_WORLD);  //j manda ad other_rank il miglior percorso
                    MPI_Recv(mesg2.data(), mesg2.size(), MPI_INTEGER, other_rank, tag2, MPI_COMM_WORLD, &stat2); //j si prepara per ricevere dal other_rank il miglior percorso
                }

                if (rank == other_rank)
                {
                    MPI_Recv(mesg1.data(), mesg1.size(), MPI_INTEGER, j, tag1, MPI_COMM_WORLD, &stat1); //*prima Recv di Send per evitare punti morti
                    MPI_Send(mesg2.data(), mesg2.size(), MPI_INTEGER, j, tag2, MPI_COMM_WORLD);
                }
                // AVVIENE LO SCAMBIO DEI MIGLIORI PERCORSI
                if (rank == j)
                {
                    generation[0] = mesg2;
                }

                if (rank == other_rank)
                {
                    generation[0] = mesg1;
                }
            }
            generation = sorting(generation, points); // perchè ho fatto degli scambi!!
        } //* la parte appena conclusa (quello della comunicazione tra cores) è stata commentata per indagare il caso di ricerca parallela ma senza comunicazione: tutto il resto, è rimasto invariato!

        ofstream WriteBestMetric;
        string filename = "core" + to_string(rank) + "_best_metric.txt";
        WriteBestMetric.open(filename, ios::app); //! appende --> NON SOVRASCRIVE I FILE PRE-ESISTENTI --> CANCELLALI PRIMA (se vuoi)!!!
        WriteBestMetric << i << "\t" << metric(generation[0], points) << endl;
        WriteBestMetric.close();

        if (rank == 0)
            cout << "Iterazione numero:\t" << i << endl;
    }

    vector<vector<int>> best_of_all(size);
    vector<int> accessorio(Npoints); // vector che conterra' il miglior percorso di ogni core

    for (int i = 0; i < size; i++)
    {
        if (rank == i)
        {
            accessorio = generation[0];
        }
        MPI_Bcast(accessorio.data(), Npoints, MPI_INTEGER, i, MPI_COMM_WORLD); // così tutti i cores (in particolar modo lo 0) sanno chi è accessorio (di volta in volta, il miglior percorso dell'i-esimo core)
        if (rank == 0) // perchè è un'operazioe che deve fare un solo core (lo 0 per convenzione)
        {
            best_of_all[i] = accessorio; // best_of_all contiene il miglior percorso di ogni core
        }
    }

    if (rank == 0) // deve farlo lo 0 (l'unico che sa cosa contiene best_of_all)
    {
        best_of_all = sorting(best_of_all, points); // la prima posizione del vettore contiene il miglior percorso di tutti i cores

        save_path(best_of_all[0], points, out); //[0] poichè voglio il cammino migliore
    }
    out.close();

    generation.clear();
    generation.shrink_to_fit();

    MPI_Finalize();

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