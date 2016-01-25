#include <iostream>
#include <stdio.h>
#include <math.h>
#include "neworkC.h"
#include "vcworkC.h"

extern"C" {
  void neutev_(int* nu_id, float e[], int* iz, int* ia);
}

extern"C" {
  void neutxs_(int* nu_in,float* e_in,
               int* z_in, int* a_in, int* mode_in,
               float* totxs, int ierr);
}

int main(){

  int nu_id = 14;
  int iz = 29;
  float  a = 64;
  int ia   = (int)a;
  int imode = 0;
  float xs = 0.0;
  int ierr;

  float e[3];
  e[0] = 0. ;
  e[1] = 0. ;
  e[2] = 20.8715964;

  float etot = e[2];

  for (int i=0; i<1000; i++){

    std::cout << i << std::endl;

    std::cout << "getting cross section " << std::endl;
    neutxs_(&nu_id, &etot, &iz, &ia, &imode, &xs);
	std::cout << xs << std::endl;

    std::cout << "getting event " << e << std::endl;
    neutev_(&nu_id, e, &iz, &ia, ierr);

    std::cout << std::endl << "Event number " << i << std::endl;
    std::cout << "Interaction mode " << nework_.modene << std::endl;
    std::cout << "Number of particles. Primary: 2 + " << nework_.numne-2 <<
      "  Secondary 2 + " << vcwork_.nvc-2 << std::endl;

    std::cout << "From common block nework (momenta in GeV/c): " << std::endl;
    for (int j=0; j<nework_.numne; j++){

      std::cout << "  Flags for particle ipne/iorgene/icrnne/iflgne " << j << " :\t " << nework_.ipne[j] <<
	"\t " << nework_.iorgne[j] << "\t " << nework_.icrnne[j] <<
	"\t " << nework_.iflgne[j] << std::endl;
      //      std::cout << "  Momentum x, y, z : " << nework_.pne[j][0] << " " <<
      // nework_.pne[j][1] << " " << nework_.pne[j][2] << std::endl;

    }

    std::cout << "From common block vcwork (momenta in MeV/c): " << std::endl;

    for (int j=0; j<vcwork_.nvc; j++){

      std::cout << "  Flags for particle: ipvc/iorgvc/icrnvc/iflgvc " << j << " :\t " << vcwork_.ipvc[j] <<
	"\t " << vcwork_.iorgvc[j] << "\t " << vcwork_.icrnvc[j] <<
	"\t " << vcwork_.iflgvc[j] << std::endl;
      std::cout << "  Momentum x, y, z : " << vcwork_.pvc[j][0] << " " <<
	vcwork_.pvc[j][1] << " " << vcwork_.pvc[j][2] << std::endl;

      std::cout << "  Position x, y, z, t : " << vcwork_.posivc[j][0] << " " <<
	vcwork_.posivc[j][1] << " " << vcwork_.posivc[j][2] << "\t" << vcwork_.timvc[j] << std::endl;

    }

    std::cout << std::endl;
    std::cout << "Enter dummy integer to continue: ";

    int dummy;
    std::cin >> dummy;

    double pneu = sqrt(vcwork_.pvc[1][0]*vcwork_.pvc[1][0] +
		      vcwork_.pvc[1][1]*vcwork_.pvc[1][1] +
		       vcwork_.pvc[1][2]*vcwork_.pvc[1][2] ); 
    int i;
    if(pneu >  225){
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "High momentum nucleon " << std::endl;
      std::cout << vcwork_.pvc[1][0] << "\t"
		<< vcwork_.pvc[1][1] << "\t"
		<< vcwork_.pvc[1][2] << "\t"
		<< pneu
		<< std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
      std::cin >> i;
    }
  }

}

