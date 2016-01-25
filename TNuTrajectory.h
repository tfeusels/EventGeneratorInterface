//////////////////////////////////////////////////////////////////////////
// 
// TNuTrajectory.h
// Main algorithm for determining interaction vertex from neutrino vector 
// information and detector geometry.
//
// The object is initalized by setting a TGeoManager, a top volume and a 
// starting z. It assumed that the neutrinos emerge from a plane perpendicular
// to the z-direction (as they are reported from jnubeam)
//
// The Swim method determines for a given neutrino vector the pathlength
// through each volume. By calculating the neutrino cross section based on 
// the material in each volume, the interaction probability in each volume
// is calculated.
//
// Based on this information, one can decide whether an interaction 
// actually occurs. This is determined by the DrawInteraction method
// based on the rejection method. 
// Where and on what kind of target nuclues it occurs handled 
// by the DrawVertex method.
// 
// To maximize the efficiency of the generation, the user should set
// fMaxIntProb, the maximum interaction probability of any neutrino vector
// for all the neutrino vectors in use. This can be done by running the
// event_rate executable. 
//////////////////////////////////////////////////////////////////////////

#include "TObjArray.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"
#include "TGeoNode.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TNuVector.h"
#include "TNuVertex.h"
#include "SystemOfUnits.h"
#include  <iostream>
#include  <vector>

//Generator specific Total Xsec headers
#include "TNeutXsec.h"
#include "TGiBUUXsec.h"

class TNuTrajectory
{

public:

  // constructors
  TNuTrajectory();
  
  // destructor
  virtual ~TNuTrajectory();
  
  // methods: 
  void Clear();
  void ClearPathInfo();

  // set random seed
  void SetRandomSeed(Int_t seed);

  // set geometry manager
  void SetGeoManager(TGeoManager* theGeoManger);
  void SetGlobalVol (TGeoVolume*  theGlobalVol);
  void SetStartZ(Double_t z);

  void SetMaxIntProb(Double_t d);
  void SetBeamWeight(Double_t d);

  // record trajectory information
  void AddPathMaterial(TGeoMaterial* theMaterial); // just for this path
  void AddMaterial(TGeoMaterial* theMaterial);     // all materials enountered ever.
  void AddNode(TGeoNode* theNode);
  void AddPathLength(Double_t thePathlength);

  // set up  scattering coefficient calculation: 
  //void     SetXsecFunc(TNeutXsec* theXsecFunc);
  //void     SetXsecFunc(TXsec* theXsecFunc);
  void     SetXsecFunc(std::string nuGenerator);
  void     ExtractNuclei();

  void     SetNuType(Int_t i);
  void     SetNuEnergy(Double_t e);

  void     CalcXsecs();
  void     CalcScatCoefficients();

  // swim the neutrino through the geometry
  void     TranslateToLocal(TGeoNode* theNode, TGeoNode* theVolume);
  void     Swim(TNuVector* theNuVector);

  // retrieve material information 
  TGeoMaterial* GetMaterial(Int_t i);
  TGeoElement*  GetNucleus(Int_t i);

  // retrieve trajectory information 
  TGeoNode*     GetNode(Int_t i);
  Double_t      GetPathLength(Int_t i);
  TGeoMaterial* GetPathMaterial(Int_t i);
  Double_t      GetPathIntProb(Int_t i);

  Double_t      GetXsec(TGeoElement* theNucleus);
  Double_t      GetScatCoefficient(TGeoMaterial* theMaterial);

  Double_t      GetMaxIntProb();
  Double_t      GetSumIntProb();
  void          CalcSumIntProb();

  Int_t         GetNMaterials();
  Int_t         GetNNuclei();
  Int_t         GetNNodes();

  Double_t       GetRawProbability();
  Bool_t         DrawInteraction();
  TLorentzVector DrawVertexFromPath(TGeoElement*& theTarget);
  TNuVertex*     DrawVertexFromPath();

protected:

  TGeoManager             *fGeoManager;      // Geometry manager  
  TGeoVolume              *fGlobalVol;       // the global volume
  Double_t                 fStartZ;          // displacement 
  TObjArray               *fPathMatList;     // list of materials for this path
  TObjArray               *fNodeList;        // list of nodes
  std::vector<Double_t>    fPathLengths;     // list of Pathlengths
  std::vector<Double_t>    fPathIntProb;     // list of int prob for each path 

  TObjArray               *fNucList;    // list of nuclei


  std::vector<Double_t>    fXsecs;      // list of neutrino cross sections
                                        // indexed by nuclei

  TObjArray               *fMatList;    // list of materials ever encountered
  std::vector<Double_t>    fScatCoeffs; // list of scattering coefficients
                                        // indexed by materials


  TXsec*               fXsecFunc;   // Function that calculates the 
                                        // neutrino cross section as a func
                                        // of E, nuType, target nucleus

  Int_t    fNuType;
  Double_t fNuEnergy;
  Double_t fMaxIntProb;

  TRandom3 fRandom;
  TVector3 fStartPoint;
  TVector3 fDir;
  Double_t fSumIntProb;
  Double_t fBeamWeight;
  
};
