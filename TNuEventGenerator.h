#ifndef TNUEVENTGENERATOR_H
#define TNUEVENTGENERATOR_H

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

class TNuEventGenerator {

public: 

  // "NE" accessors and modifiers correspond to the particle stack
  // prior to the effect of secondary interactions in the target nucleus
  // i.e. the outgoing particles of the "bare" neutrino/nucleon interactions.

  // "VC" accessors and modifiers correspond to the particle stack
  // AFTER the effect of secondary interactions in the target nucleus.
 
  TNuEventGenerator(){
    fPNE = new float[3];
    fPVC = new float[3];
    fPosVC = new float[4];
    fPosVert = new float[3];
    fDir = new float[3];
  };

  ~TNuEventGenerator(){
    delete fPNE;
    delete fPVC;
    delete fPosVC;
    delete fPosVert;
    delete fDir;
    fPNE = 0;
    fPVC = 0;
    fPosVC = 0;
    fPosVert = 0;
    fDir = 0;
  };

  // neutrino event generator information:
  // neutrino type, momentum vector in GeV/c, Z and A of target nucleus
  virtual void getVector(int nutype, float* pgevc, int Z, int iA, int &ierr) = 0;

  virtual int getInteractionMode() = 0;
  virtual int getNNEParticles() = 0;
  virtual int getNVCParticles() = 0;


  // Neutrino Interaction vertex in nucleus in fm relative to center 
  // of nucleus. 
  virtual float *getVCPos(int ip) = 0;


  virtual float* getNEMomentum(int ip) = 0;
  virtual float* getVCMomentum(int ip) = 0;
  virtual float getSecMass(int ip) = 0;

  virtual int getNEID(int ip) = 0;
  virtual int getNEOrigin(int ip) = 0;
  virtual int getNEChase(int ip) = 0;
  virtual int getNEFlag(int ip) = 0;

  virtual int getVCID(int ip) = 0;
  virtual int getVCOrigin(int ip) = 0;

  virtual int getVCChase(int ip) = 0;
  virtual int getVCFlags(int ip) = 0;


  // FSI history information
  virtual int getNFSVertices() = 0;
  virtual int getFSIntFlag(int iv) = 0;
  virtual float*  getFSPos(int iv) = 0;

  virtual int getNFSParticles() = 0;
  virtual int getFSID(int ip) = 0;
  virtual int getFSVertexInit(int ip) = 0;
  virtual int getFSVertexFin(int ip) = 0;
  virtual float getFSMomentumLab(int ip) = 0;
  virtual float getFSMomentumNuc(int ip) = 0;
  virtual float*  getFSDir(int ip) = 0;


  // Cross section and bound proton info
  virtual float getTotcrs() = 0;
  virtual int   getIbound() = 0;

  // Coherent Cross section info
  virtual float getCrsX() = 0;
  virtual float getCrsY() = 0;
  virtual float getCrsZ() = 0;
  virtual float getCrsPhi() = 0;

protected:
  
  float* fPNE;
  float* fPVC;
  float* fPosVC;

  float* fPosVert;
  float* fDir;
   
  int   fZ, fIA;

};


#endif
