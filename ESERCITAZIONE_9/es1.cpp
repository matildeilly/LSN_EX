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
#include <algorithm>

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
{ // controlla che il vettore sia un cammino accettabile
    for (int i = 0; i < V.size(); i++)
    {
        int count = 0;
        for (int j = 0; j < V.size(); j++)
        {
            if (j == 0 && V[j] != 0) // controllo che il primo elemento sia sempre 0
            {
                return false;
            }
            if (V[j] == i)
            {
                count++;
            }
        }
        if (count != 1) // controllo che ogni percorso contenga ogni città solo una volta
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
            checked_gen.push_back(gen[i]);  //se il percorso va bene lo aggiungo alla nuova generazione
        assert(check(gen[i]) && "Invalid path generated after check");
    }
    if (checked_gen.size() % 2 != 0)
    {                                                            // voglio che la dim sia pari (vedi crossing...)
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
    return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2)); //è la consueta distanza tra due punti
}

double metric(vector<int> A /*percorso di cui calcolo la lunghezza*/, vector<position> points /*vettore che consente di assegnare ad ogni città un numero (dato dall'indice della città in points)*/) //*calcola la lunghezza di un percorso 
{
    int dim = A.size();
    vector<position> V(dim);
    for (int i = 0; i < dim; i++)
    {
        V[i] = points[A[i]]; //V è la stessa cosa di A MA ho convertito i numeri in position
    }
    double add = 0.; 
    for (int i = 0; i < (dim - 1); i++)
    { //(dim - 1) così all'ultimo step ho la distanza tra il penultimo e l'ultimo elemento
        add += dist(V[i], V[i + 1]);
    }
    add += dist(V[0], V[dim - 1]); // distanza tra l'ultimo elemento e il primo (devo ritornare al punto di partenza!)
    return add;
}

int identifica_pos_in_vec(double val, vector<double> V) //restituisce il primo indice del vettore V in cui compare il valore val
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
        new_vec[i] = V[identifica_pos_in_vec(metriche_to_sort[i], metriche)];   //new_vec è la generazione ordinata (contiene gli stessi elementi di V ma ordinato per metriche crescenti) 
    }
    return new_vec;
}

vector<int> select(vector<vector<int>> V, Random *rnd)  //seleziono con maggior probabilità i primi percorsi (i più brevi)
{
    int M = V.size();
    double p = 2;
    int j = (int)(M * pow(rnd->Rannyu(), p));
    return V[j];
}

void mutation_1(vector<int> &V, Random *rnd)    //scambio tra loro due città casualmente
{
    int dim = V.size();
    int index_1 = (int)rnd->Rannyu(1., dim);
    int index_2 = (int)rnd->Rannyu(1., dim);
    swap(V[index_1], V[index_2]);
    assert(check(V) && "Invalid path generated after mutation_1");
}

void mutation_2(vector<int> &V, Random *rnd) // shifto n città a partire dalla m-esima fino in fondo
{
    int N = V.size();
    int m = (int)rnd->Rannyu(1, N - 1); // indice da cui partire
    int n = (int)rnd->Rannyu(1, N - m); // numero di città da prendere

    vector<int> temp(n); // temp è un vettore accessorio...

    for (int i = 0; i < n; ++i)
    { // in temp salvo il blocco che voglio spostare
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

void mutation_3(vector<int> &V, Random *rnd) // swappo due blocchi
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

void mutation_4(vector<int> &V, Random *rnd)    // scambio l'ordine delle città in un dato blocco 
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

void kaos(vector<vector<int>> &gen, Random *rnd, double prob)   // faccio una mutazione (scelta a caso) con probabilità prob
{ // prob è la probabilità di FARE la mutazione 
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

vector<vector<int>> initialize_generation(Random *rnd, int dim_gen, int number_of_points, int number_of_permut) // la prima generazione la ottengo a partire da un vettore ordinato di interi per poi scambiare tra loro due città più volte a caso (tramite mutation_1)
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

bool search(int value, const vector<int> &V) //ci dice se all'interno di V c'è value
{ // const per assicurarci che search non modifichi V e & (by reference) per efficienza
    for (int i = 0; i < V.size(); i++)
    {
        if (V[i] == value)
            return true;
    }
    return false;
}

vector<int> update_small(vector<int> &small_parent, const vector<int> &other_parent) // aggiunge tutti gli elementi di other_parent che non sono presenti in small_parent nell'ordine in cui compaiono in other_parent
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

int identifica_pos_in_gen(vector<int> V, vector<vector<int>> gen) //restituisce il primo indice del vettore gen in cui compare il vettore V
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
{ // prob di FARE il crossing 
    int dim_gen = gen.size();
    vector<vector<int>> next_gen;
    //voglio che la generazione successiva abbia lo stesso numero di elementi di quella corrente (gen) 
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

        if (discriminator > prob) //in questo caso non avviene il crossing --> i genitori vanno nella nuova generazione
        {
            next_gen.push_back(parent_1);
            next_gen.push_back(parent_2);
        }
        else // faccio il crossing
        {
            int N = parent_1.size();              // per forza la stessa di parent_2...
            int cut_off = (int)rnd->Rannyu(1, N); // assicurarsi che il taglio non sia all'inizio o alla fine
            vector<int> small_parent_1(parent_1.begin(), parent_1.begin() + cut_off); // ha gli stessi elementi dei genitori fino al taglio
            vector<int> small_parent_2(parent_2.begin(), parent_2.begin() + cut_off); // ha gli stessi elementi dei genitori fino al taglio

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

vector<position> points_circle(Random *rnd, double r, int number_of_points) // genero punti su una circonferenza di raggio r
{
    vector<position> points(number_of_points);
    for (int i = 0; i < number_of_points; i++)
    {
        double theta = rnd->Rannyu(0., 2 * M_PI);
        points[i].x = r * cos(theta);
        points[i].y = r * sin(theta);
    }
    return points;
}

vector<position> points_square(Random *rnd, double l, int number_of_points) // genero punti in un quadrato di lato 2l
{
    vector<position> points(number_of_points);
    for (int i = 0; i < number_of_points; i++)
    {
        points[i].x = rnd->Rannyu(-l, l);
        points[i].y = rnd->Rannyu(-l, l);
    }
    return points;
}

void save_path(vector<int> best_path, vector<position> points, ofstream &out)
{
    for (int i = 0; i < best_path.size(); i++)
    {
        int index = best_path[i];
        out << points[index].x << "\t" << points[index].y << endl; // salvo le città nell'ordine in cui sono in best_path
    }
    out << points[0].x << "\t" << points[0].y << endl; // devo tornare al punto di partenza
}

void save_length(vector<int> best_path, vector<position> points, ofstream &out, int index) // savlo la lunghezza del miglior percorso
{
    out << index << "\t" << metric(best_path, points) << endl;
}

void save_mean_length(vector<vector<int>> gen, vector<position> points, ofstream &out, int index) // salvo la media delle lunghezze per la prima metà (quella dei migliori) dei percorsi nella generazione
{
    double sum = 0.;
    for (int i = 0; i < gen.size() / 2.; i++)
    {
        sum += metric(gen[i], points);
    }
    out << index << "\t" << sum / (gen.size() / 2.) << endl;
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

    ofstream ciout("circle_path.txt");
    ofstream sqout("square_path.txt");
    ofstream cmout("circle_best_metric.txt");
    ofstream smout("square_best_metric.txt");
    ofstream cmmout("circle_mean_metric.txt");
    ofstream smmout("square_mean_metric.txt");

    int Npoints = 34;
    double r = 1;   //raggio cerchio
    double l = 0.5; //metà del lato del quadrato
    vector<position> circle = points_circle(&rnd, r, Npoints);
    vector<position> square = points_square(&rnd, l, Npoints);
    int dim_gen = 150; //dimensione della generazione
    vector<vector<int>> generation = initialize_generation(&rnd, dim_gen, Npoints, 100);
    //* generation = check_gen(generation);

    int Niter = 500; //numero di iterazioni
    double prob_cross = 0.7; //probabilità di crossover
    double prob_mut = 0.2; //probabilità di mutazione
    vector<vector<int>> generation_circle = sorting(generation, circle);
    vector<vector<int>> generation_square = sorting(generation, square);
    for (int i = 0; i < Niter; i++)
    {
        generation_circle = crossing(generation_circle, &rnd, prob_cross);
        kaos(generation_circle, &rnd, prob_mut); //faccio le mutazioni
        //* generation_circle = check_gen(generation_circle);
        generation_circle = sorting(generation_circle, circle);

        generation_square = crossing(generation_square, &rnd, prob_cross);
        kaos(generation_square, &rnd, prob_mut);
        //* generation_square = check_gen(generation_square);
        generation_square = sorting(generation_square, square);

        save_length(generation_circle[0], circle, cmout, i); //[0] poichè voglio il cammino migliore
        save_length(generation_square[0], square, smout, i); //[0] poichè voglio il cammino migliore
        save_mean_length(generation_circle, circle, cmmout, i);
        save_mean_length(generation_square, square, smmout, i);
    }
    save_path(generation_circle[0], circle, ciout); //[0] poichè voglio il cammino migliore
    save_path(generation_square[0], square, sqout); //[0] poichè voglio il cammino migliore

    ciout.close();
    sqout.close();
    cmout.close();
    smout.close();
    cmmout.close();
    smmout.close();

    generation_circle.clear();
    generation_circle.shrink_to_fit();
    generation_square.clear();
    generation_square.shrink_to_fit();

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
