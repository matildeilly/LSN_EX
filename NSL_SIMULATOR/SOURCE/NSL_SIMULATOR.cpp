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
#include "system.h"
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{

  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);

  //  ofstream out("instant_en_solid_out.txt");

  for (int i = 0; i < 10000 /* (SOLID) 30000 (LIQUID)250000(GAS)*/; i++)
  {
    SYS.step();
  } // * * EQUILIBRAZIONE ESERCITAZIONE 4

  /*for (int i=0; i<200; i++){
      SYS.step();
    }  // * * EQUILIBRAZIONE ESERCITAZIONE 6 */

  //  for (int i=0; i<1000/*every case*/; i++){
  //    SYS.step();
  //  } // * * EQUILIBRAZIONE ESERCITAZIONE 7
  

  // * * EQUILIBRIUM: è stata usata in tutte le esercitazioni per monitorare la grandezza di interesse al fine di capire dopo quanti steps la simulazione è in equilibrio. Una volta capito ciò, l'equilibrazione è stata performata tramite degli steps (vedi sopra)

  /*  for (int i=0; i<500000; i++){                 
      out<<i<<"\t"<<SYS.equilibrium()<<endl; }*/  

  for (int i = 0; i < SYS.get_nbl(); i++)
  { // loop over blocks
    for (int j = 0; j < SYS.get_nsteps(); j++)
    { // loop over steps in a block
      SYS.step();
      SYS.measure();
      if (j % 10 == 0)
      {
        //        SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
        nconf++;
      }
    }
    SYS.averages(i + 1);
    SYS.block_reset(i + 1);
  }
  SYS.finalize();

  //  out.close();
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
