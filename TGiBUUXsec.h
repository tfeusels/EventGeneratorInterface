/************************************************************
 * TNeutXsec.h
 *
 * GiBUU interface of total cross section calculation.
 * The access to the total cross section is done through
 * a wrapper around the fortran function neutxs(.F) linked
 * in the GNUMakefile to event_rate.o
 *
 ************************************************************/

#include "SystemOfUnits.h"
#include "TXsec.h"

extern "C" {
  void gibuuxs_(int* nu_in,float* e_in,
               int* z_in, int* a_in, int* mode,
               float* totxs);
}

class TGiBUUXsec : public TXsec {

public:

 TGiBUUXsec():TXsec() { 
    fExponentCm2 = -38;
  };

 ~TGiBUUXsec(){ };

  float getTotalXsec(int nutype, float e, int z, int a) 
  {  
    float xsec;
    float egev = e/CLHEP::GeV;
    int   ztemp = z;
    int   atemp = a;
    int   mode  = 0; //finalstate_ID, 0 = total!
                     //TODO: NC vs CC (how is it done in NEUT??)

    gibuuxs_(&nutype, &egev, &ztemp, &atemp, &mode, &xsec);

    return xsec;
  };


};
