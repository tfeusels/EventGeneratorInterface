#ifndef TNEUTEVENTGENERATOR_H
#define TNEUTEVENTGENERATOR_H

//////////////////////////////////////////////////////////////////////////////
// TNeutEventGenerator.h
//
// Calls the neutev function to generate neutrino event vectors 
// Provides access to the common blocks (nework, vcwork, fsihist) which specify
// the final state of the interactions
//
// nework: Particle state particles existing prior to final state interactions
//         nework quantities are accessed via "NExxx" methods
//
// vcwork: Particle state including propagation through target nucleus
//         vcwork quantities are accessed via "VCxxx" methods
//
// fsihist: Hadron (currently only pion) interaction history through nucleus
//          fsihist quantities are accessed via "FSxxx" methods
//
// necrs: Stores calculated total (future: differential) cross sections
//
// posinnuc: Stores if target proton is bound or not. And position within 
//           nucleus (not yet implemented here).
//
// "Native" neut units are preserved without conversion to CLHEP
//////////////////////////////////////////////////////////////////////////////

// common block structures from Neut
#include "vcworkC.h"
#include "neworkC.h"
#include "fsihistC.h"
#include "neutcrsC.h"
#include "posinnucC.h"
#include "TNuEventGenerator.h"

extern"C" {
  void neutev_(int* nu_id, float e[], int* iz, int* ia, int* ierr);
}

class TNeutEventGenerator: public TNuEventGenerator{

public: 

  // "NE" accessors and modifiers correspond to the particle stack
  // prior to the effect of secondary interactions in the target nucleus
  // i.e. the outgoing particles of the "bare" neutrino/nucleon interactions.

  // "VC" accessors and modifiers correspond to the particle stack
  // AFTER the effect of secondary interactions in the target nucleus.
 
 TNeutEventGenerator(): TNuEventGenerator(){
    //anything Neut specific?
  };

  ~TNeutEventGenerator(){
    //anything Neut specific?
  };

  // neutrino event generator information:
  // neutrino type, momentum vector in GeV/c, Z and A of target nucleus
  void getVector(int   nutype, float* pgevc,
		 int   Z,      int    iA, int &ierr){
    fZ  =  Z;
    fIA = iA;
    neutev_(&nutype, pgevc, &fZ, &fIA, &ierr);
  };

  int getInteractionMode() { return nework_.modene;}
  int getNNEParticles()   { return nework_.numne; }
  int getNVCParticles()   { return vcwork_.nvc;   }


  // Neutrino Interaction vertex in nucleus in fm relative to center 
  // of nucleus. 
  float *getVCPos(int ip) {
    fPosVC[0] = vcwork_.posvc[0];
    fPosVC[1] = vcwork_.posvc[1];
    fPosVC[2] = vcwork_.posvc[2];
    fPosVC[3] = vcwork_.timvc[ip]; // may have no meaning currently
    return fPosVC;
  }


  float* getNEMomentum(int ip) {
    // momentum returned in GeV/c
    // convert to CLHEP stanard MeV/c
    fPNE[0] = 1000.0 * nework_.pne[ip][0];
    fPNE[1] = 1000.0 * nework_.pne[ip][1];
    fPNE[2] = 1000.0 * nework_.pne[ip][2];
    return fPNE;
  };

  float*  getVCMomentum(int ip) {
    // momentum returned in MeV/c from the VCWORK common block
    // this is CLHEP default units so don't convert
    fPVC[0] = vcwork_.pvc[ip][0];
    fPVC[1] = vcwork_.pvc[ip][1]; 
    fPVC[2] = vcwork_.pvc[ip][2];
    return fPVC;
  };

  float  getSecMass(int ip){
    // mass is returned in MeV/c^2 from VCWORK block
    // this CLHEP default so don't need to convert.
    return vcwork_.amasvc[ip];
  }


  int getNEID(int ip)     { return nework_.ipne[ip];   }
  int getNEOrigin(int ip) { return nework_.iorgne[ip]; }
  int getNEChase(int ip)  { return nework_.icrnne[ip]; }
  int getNEFlag(int ip)   { return nework_.iflgne[ip]; }

  int getVCID(int ip)     { return vcwork_.ipvc[ip];   }

  // iorgvc follows Fortran numbering convention starting at 1.
  // This should be fine when refering to the actual C++ index 
  // since we are adding the nucleus into slot 1.
  int getVCOrigin(int ip) { return vcwork_.iorgvc[ip]; }

  int getVCChase(int ip)  { return vcwork_.icrnvc[ip]; }
  int getVCFlags(int ip)  { return vcwork_.iflgvc[ip]; }


  // FSI history information
  int getNFSVertices()   { return fsihist_.nvert;   }
  int getFSIntFlag(int iv)   { return fsihist_.iflgvert[iv];   }
  float*  getFSPos(int iv) {
    fPosVert[0] = fsihist_.posvert[iv][0];
    fPosVert[1] = fsihist_.posvert[iv][1]; 
    fPosVert[2] = fsihist_.posvert[iv][2];
    return fPosVert;
  };

  int getNFSParticles()   { return fsihist_.nvcvert;   }
  int getFSID(int ip)     { return fsihist_.ipvert[ip];   }
  int getFSVertexInit(int ip)     { return fsihist_.iverti[ip];   }
  int getFSVertexFin(int ip)     { return fsihist_.ivertf[ip];   }
  float getFSMomentumLab(int ip) { return fsihist_.abspvert[ip]; }
  float getFSMomentumNuc(int ip) { return fsihist_.abstpvert[ip]; }
  float*  getFSDir(int ip) {
    fDir[0] = fsihist_.dirvert[ip][0];
    fDir[1] = fsihist_.dirvert[ip][1]; 
    fDir[2] = fsihist_.dirvert[ip][2];
    return fDir;
  };


  // Cross section and bound proton info
  float getTotcrs() { return neutcrscom_.totcrsne; }
  int   getIbound() { return posinnuc_.ibound; }

  // Coherent Cross section info
  float getCrsX() { return neutcrscom_.crsx; }
  float getCrsY() { return neutcrscom_.crsy; }
  float getCrsZ() { return neutcrscom_.crsz; }
  float getCrsPhi() { return neutcrscom_.crsphi; }

};


#endif
