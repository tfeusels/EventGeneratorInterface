#ifndef TNuVector_h
#define TNuVector_h

#include "TVector3.h"
#include "SystemOfUnits.h"

class TNuVector {

public:

  TNuVector()  {fNumIntSteps=0;};
  ~TNuVector() {};


  void      setNuFileName(TString s){fNuFileName = s;}
  TString   getNuFileName(){ return fNuFileName; }

  void      setEntry(Int_t i) 
            { fEntry = i; }
  Int_t     getEntry()
            { return fEntry; }

  void      setNuFluxEntry(Int_t i) 
            { fNuFluxEntry = i; }
  Int_t     getNuFluxEntry()
            { return fNuFluxEntry; }

  void      setNuEnergy(Double_t e)    
            { fNuEnergy = e;    }
  Double_t  getNuEnergy()              
            { return fNuEnergy/CLHEP::MeV; }

  void      setDecayParticleID(Int_t i)        
            { fDecayParticleID = i; }
  Int_t     getDecayParticleID()               
            { return fDecayParticleID; }

  void      setDecayMode(Int_t i)              
            { fDecayMode = i; }
  Int_t     getDecayMode()                     
            { return fDecayMode; }

  void      setNuType(Int_t i)         
            { fNuType = i;      }
  Int_t     getNuType()                
            { return fNuType;   }

  void      setDecayPos(Double_t x, Double_t y,Double_t z)  
            { fDecayPos.SetXYZ(x,y,z); }
  TVector3  getDecayPos()                  
            { return fDecayPos;        }

  void      setDecayMomentum(Double_t x, Double_t y,Double_t z) 
            { fDecayMomentum.SetXYZ(x,y,z); }
  TVector3  getDecayMomentum()                 
            { return fDecayMomentum;        }

  void      setCosToBeamAtDecay(Float_t f)
            { fCosToBeamAtDecay = f;    }
  Float_t   getCosToBeamAtDecay()
            { return fCosToBeamAtDecay; }


  void      setNuWeight(Float_t f){ fNuWeight = f;    }
  Float_t   getNuWeight()         { return fNuWeight; }


  void      setNvtx0(Int_t i){ fNvtx0 = i;    }
  Int_t     getNvtx0()         { return fNvtx0; }




  void      setProductionPos(Double_t x, Double_t y,Double_t z) 
            { fProductionPos.SetXYZ(x,y,z); }
  TVector3  getProductionPos()                 
            { return fProductionPos;        }

  void      setProductionMomentum(Double_t x, Double_t y,Double_t z) 
            { fProductionMomentum.SetXYZ(x,y,z); }
  TVector3  getProductionMomentum()                
            { return fProductionMomentum;        }

  void      setCosToBeamAtProduction(Float_t f)
            { fCosToBeamAtProduction = f;    }
  Float_t   getCosToBeamAtProduction()
            { return fCosToBeamAtProduction; }


  void      setNuRnu(Double_t r) 
            { fNuRnu = r; }
  Double_t  getNuRnu()                 
            { return fNuRnu;        }

  void      setNuPos(Double_t x, Double_t y,Double_t z) 
            { fNuPos.SetXYZ(x,y,z); }
  TVector3  getNuPos()                 
            { return fNuPos;        }


  void      setNuDir(Double_t x, Double_t y,Double_t z) 
            { fNuDir.SetXYZ(x,y,z); }

  TVector3  getNuDir()                 
            { return fNuDir;        }


  void      setIdfd(Int_t i){ fIdfd = i;    }
  Int_t     getIdfd()         { return fIdfd; }



  void      setPrimaryParticleID(Int_t i)
            { fPrimaryParticleID = i; }
  Int_t     getPrimaryParticleID()
            { return fPrimaryParticleID; }

  void      setPrimaryParticlePos(Double_t x, Double_t y, Double_t z)
            { fPrimaryParticlePos.SetXYZ(x,y,z); }
  TVector3  getPrimaryParticlePos()
            { return fPrimaryParticlePos; }

  void      setPrimaryParticleMomentum(Double_t x, Double_t y, Double_t z)
            { fPrimaryParticleMomentum.SetXYZ(x,y,z); }
  TVector3  getPrimaryParticleMomentum()
            { return fPrimaryParticleMomentum; }



  // 2011a
  void      setNumIntSteps(Int_t i) 
            { fNumIntSteps = i; }
  Int_t     getNumIntSteps()                 
            { return fNumIntSteps;        }

  void      setIntermediateParticleMomentum(Float_t x[], Float_t y[], Float_t z[]) {
              for (Int_t i=0; i<fNumIntSteps; i++)
		fIntermediateParticleMomentum[i].SetXYZ(x[i]*CLHEP::GeV,y[i]*CLHEP::GeV,z[i]*CLHEP::GeV); 
  }
  TVector3  getIntermediateParticleMomentum(Int_t i)
            { return fIntermediateParticleMomentum[i]; }


  void      setCosToBeamOfIntermediate(Float_t f[]) {
              for (Int_t i=0; i<fNumIntSteps; i++)
		fCosToBeamOfIntermediate[i] = f[i]; 
  }
  Float_t   getCosToBeamOfIntermediate(Int_t i)
            { return fCosToBeamOfIntermediate[i]; }


  void      setIntermediateParticlePos(Float_t x[], Float_t y[], Float_t z[]) {
              for (Int_t i=0; i<fNumIntSteps; i++)
		fIntermediateParticlePos[i].SetXYZ(x[i]*CLHEP::cm,y[i]*CLHEP::cm,z[i]*CLHEP::cm); 
  }
  TVector3  getIntermediateParticlePos(Int_t i)
            { return fIntermediateParticlePos[i]; }


  void      setIntermediateParticleID(Int_t id[]) {
              for (Int_t i=0; i<fNumIntSteps; i++)
		fIntermediateParticleID[i] = id[i]; 
  }
  Int_t     getIntermediateParticleID(Int_t i)
            { return fIntermediateParticleID[i]; }

  void      setIntermediateParticleMec(Int_t mec[]) {
              for (Int_t i=0; i<fNumIntSteps; i++)
		fIntermediateParticleMec[i] = mec[i]; 
  }
  Int_t     getIntermediateParticleMec(Int_t i)
            { return fIntermediateParticleMec[i]; }


  void      setNuEnergySK(Double_t e)    
            { fNuEnergySK = e;    }
  Double_t  getNuEnergySK()              
            { return fNuEnergySK; }

  void      setNuWeightSK(Float_t f){ fNuWeightSK = f;    }
  Float_t   getNuWeightSK()         { return fNuWeightSK; }

  void      setNuWeightAccept(Float_t f){ fNuWeightAccept = f;    }
  Float_t   getNuWeightAccept()         { return fNuWeightAccept; }


  void      setNuGmat(Int_t gmat[]) {
              for (Int_t i=0; i<fNumIntSteps; i++)
		fNuGmat[i] = gmat[i]; 
  }
  Int_t     getNuGmat(Int_t i)
            { return fNuGmat[i]; }

  void      setNuGdistc(Float_t f[]) {
              for (Int_t i=0; i<fNumIntSteps; i++)
		fNuGdistc[i] = f[i]*CLHEP::cm; 
  }
  Float_t   getNuGdistc(Int_t i)
            { return fNuGdistc[i]; }

  void      setNuGdistal(Float_t f[]) {
              for (Int_t i=0; i<fNumIntSteps; i++)
		fNuGdistal[i] = f[i]*CLHEP::cm; 
  }
  Float_t   getNuGdistal(Int_t i)
            { return fNuGdistal[i]; }

  void      setNuGdistti(Float_t f[]) {
              for (Int_t i=0; i<fNumIntSteps; i++)
		fNuGdistti[i] = f[i]*CLHEP::cm; 
  }
  Float_t   getNuGdistti(Int_t i)
            { return fNuGdistti[i]; }

  void      setNuGdistfe(Float_t f[]) {
              for (Int_t i=0; i<fNumIntSteps; i++)
		fNuGdistfe[i] = f[i]*CLHEP::cm; 
  }
  Float_t   getNuGdistfe(Int_t i)
            { return fNuGdistfe[i]; }



  void      setNuVersion(Float_t f){fNuVersion = f;}
  Float_t   getNuVersion(){ return fNuVersion; }

  void      setNuTuneid(Int_t i){fNuTuneid = i;}
  Int_t     getNuTuneid(){ return fNuTuneid; }

  void      setNuNtrig(Int_t i){fNuNtrig = i;}
  Int_t     getNuNtrig(){ return fNuNtrig; }

  void      setNuPint(Int_t i){fNuPint = i;}
  Int_t     getNuPint(){ return fNuPint; }

  void      setNuBpos(Float_t f[]) { for (int i=0; i<2; i++) fNuBpos[i] = f[i]; }
  Float_t   getNuBpos(Int_t i) { return fNuBpos[i]; }

  void      setNuBtilt(Float_t f[]) { for (int i=0; i<2; i++) fNuBtilt[i] = f[i]; }
  Float_t   getNuBtilt(Int_t i) { return fNuBtilt[i]; }

  void      setNuBrms(Float_t f[]) { for (int i=0; i<2; i++) fNuBrms[i] = f[i]; }
  Float_t   getNuBrms(Int_t i) { return fNuBrms[i]; }

  void      setNuEmit(Float_t f[]) { for (int i=0; i<2; i++) fNuEmit[i] = f[i]; }
  Float_t   getNuEmit(Int_t i) { return fNuEmit[i]; }

  void      setNuAlpha(Float_t f[]) { for (int i=0; i<2; i++) fNuAlpha[i] = f[i]; }
  Float_t   getNuAlpha(Int_t i) { return fNuAlpha[i]; }

  void      setNuHcur(Float_t f[]) { for (int i=0; i<3; i++) fNuHcur[i] = f[i]; }
  Float_t   getNuHcur(Int_t i) { return fNuHcur[i]; }

  void      setNuRand(Int_t i){fNuRand = i;}
  Int_t     getNuRand(){ return fNuRand; }

  void      setNuRseed(Int_t i[]){fNuRseed[0] = i[0]; fNuRseed[1] = i[1];}
  Int_t     getNuRseed(Int_t i){ return fNuRseed[i]; }


private:

  TString fNuFileName;

  Int_t    fEntry, fNuFluxEntry;

  Double_t fNuEnergy;
  Int_t    fDecayParticleID, fDecayMode, fNuType;

  TVector3 fDecayPos, fDecayMomentum;
  Float_t  fCosToBeamAtDecay;
  
  Float_t  fNuWeight;
  Int_t    fNvtx0;

  TVector3 fProductionPos, fProductionMomentum;
  Float_t  fCosToBeamAtProduction;

  Double_t fNuRnu;

  TVector3 fNuPos, fNuDir;

  Int_t fIdfd;

  Int_t    fPrimaryParticleID;
  TVector3 fPrimaryParticlePos;
  TVector3 fPrimaryParticleMomentum;


  // >= 2010d
  Int_t    fNumIntSteps;
  TVector3 fIntermediateParticleMomentum[12];
  TVector3 fIntermediateParticlePos[12];
  Float_t  fCosToBeamOfIntermediate[12];
  Int_t    fIntermediateParticleID[12];
  Int_t    fIntermediateParticleMec[12]; 

  Double_t fNuEnergySK;
  Float_t  fNuWeightSK, fNuWeightAccept;

  Int_t    fNuGmat[12]; 
  Float_t         fNuGdistc[12];  
  Float_t         fNuGdistal[12]; 
  Float_t         fNuGdistti[12]; 
  Float_t         fNuGdistfe[12]; 

  // Beam parameter information (>=2010d)
  Float_t  fNuVersion;
  
  Int_t fNuTuneid;
  Int_t fNuNtrig;
  Int_t fNuPint;
  Float_t fNuBpos[2], fNuBtilt[2];
  Float_t fNuBrms[2];
  Float_t fNuEmit[2], fNuAlpha[2];
  Float_t fNuHcur[3];
  Int_t fNuRand;
  Int_t fNuRseed[2];

};
#endif
