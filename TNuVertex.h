//////////////////////////////////////////////////////////////////////////
// 
// Tnuvertex.H
// Contains Information About The Vertex Of The Generated Interaction
// Target Nucleus And 4-Position
//
//////////////////////////////////////////////////////////////////////////

#ifndef TNuVertex_h
#define TNuVertex_h

#include "TLorentzVector.h"
#include "TMath.h"
#include "TGeoElement.h"

class TNuVertex {

public:

  TNuVertex(){};
  TNuVertex(TLorentzVector theVertex,
	    TGeoElement*   theNucleus){
    setVertex(theVertex);
    setNucleus(theNucleus);
  }
  
  TNuVertex(TLorentzVector theVertex,
	    Int_t Z, Float_t A){
    setVertex(theVertex);
    setZ(Z);
    setA(A);
    setIA(TMath::Nint(A));
  }

  ~TNuVertex() {};

 
  TLorentzVector getVertex() { return fVertex; };
  Int_t    getZ()      { return fZ;    }
  Int_t    getIA()     { return fIA;   }
  Float_t  getA()      { return fA;    }
  Double_t getXsec()   { return fXsec; }

  void     setVertex(TLorentzVector vertex) { fVertex = vertex; };
  void     setVertex(Float_t t, Float_t x, Float_t y, Float_t z)
                                     { fVertex.SetXYZT(x, y, z, t); }

  void     setZ (Int_t i)              { fZ  = i; }
  void     setIA(Int_t i)              { fIA = i; }
  void     setA (Float_t f)            { fA  = f; }
  void     setNucleus(TGeoElement* theNucleus)
                                      { fZ  = theNucleus->Z(); 
					fIA = TMath::Nint(theNucleus->A()); 
					fA  = theNucleus->A();
				      }

  void     setXsec(Double_t d)        { fXsec = d; }

  int      getPDGofNucleus()          {
                                        Int_t pdgCode = atoi(Form("100%03d%03d0", fZ, fIA));
                                        return pdgCode;
                                      }

private:

  Int_t    fZ, fIA;
  Float_t  fA;
  Double_t fXsec;
  TLorentzVector fVertex;
 
};
#endif
