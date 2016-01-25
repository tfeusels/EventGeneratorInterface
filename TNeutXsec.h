#ifndef TNEUTXSEC_H
#define TNEUTXSEC_H

/************************************************************
 * TNeutXsec.h
 *
 * NEUT interface of total cross section calculation.
 * The access to the total cross section is done through
 * a wrapper around the fortran function neutxs(.F) linked
 * in the GNUMakefile to event_rate.o
 *
 ************************************************************/

#include "SystemOfUnits.h"
#include "TXsec.h"

extern "C" {
  void neutxs_(int* nu_in,float* e_in,
               int* z_in, int* a_in, int* mode,
               float* totxs);
}

class TNeutXsec : public TXsec {

public:

 TNeutXsec():TXsec() { 
    fExponentCm2 = -38;
  };

 ~TNeutXsec(){ };

  float getTotalXsec(int nutype, float e, int z, int a) 
  {  
    float xsec;
    float egev = e/CLHEP::GeV;
    int   ztemp = z;
    int   atemp = a;
    int   mode  = 0; // Choose NORMAL mode (all channels)
                     // CAN in principal through ND280Control pick any mode!!
                     // TODO : need to use this as input param then and update TXsec.h

    neutxs_(&nutype, &egev, &ztemp, &atemp, &mode, &xsec);

    return xsec;
  };

  


};
#endif
