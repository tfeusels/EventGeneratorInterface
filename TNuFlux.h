//////////////////////////////////////////////////////////////////////////
// 
// TNuFlux.h
// Accepts neutrino information from jnubeam file.
// Jnubeam files are chained together 
// TNuFlux then selects neutrinos based on:
// 1. detector location (5 for nd280, 2 for on-axis, 3/4 INGRID)
// 2. desired neutrino type (0 for all)
//
// The getVector method will read in neutrino vectors from the tree
// until an appropriate neutrino vector is found.
// The chain can be rewound if one wishes to recycle them
//////////////////////////////////////////////////////////////////////////
#ifndef TNuFlux_h
#define TNuFlux_h

#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TNuVector.h"
#include "SystemOfUnits.h"
#include <iostream>

class TNuFlux {

public:

  TNuFlux();
  ~TNuFlux();
  TNuFlux(TString fileName,TString treeKey, TString verKey, TString ndKey);

  void      initialize();
  void      addToChain(TString theFile);

  void      setTreeKey(TString theKey);
  TString   getTreeKey();

  void      setVerTreeKey(TString verKey);
  TString   getVerTreeKey();
  int       checkVerTree(TString theFile);
  Float_t   getVersion(){return version;};

  void      setNDTreeKey(TString ndKey);
  TString   getNDTreeKey();
  void      setPlanePos();

  Int_t     getNChainFluxFiles() { return nChainFluxFiles; };
  Int_t     getNHaddFluxFiles() { return nHaddFluxFiles; };

  // select specific detector
  // 0 = all, 1 = 2 km,  2 = on-axis for 2.5 deg, 3 = INGRID x, 
  // 4 = INGRID y (set to 3 for now), 5 = nd280 basket, 6 = nd280 magnet
  void     setDetector(Int_t i);
  Int_t    getDetector();

  // select neutrino type
  // 0 = all , 1 = numu, 2 = numub, 3 = nue, 4 = nueb
  void     setNuType(Int_t i);
  Int_t    getNuType();
  
  TNuVector* getVector();

  void     FillVector(TNuVector* theNuVector);
  Bool_t   ValidVector(Int_t detectorID, Int_t neutrinoID);

  Int_t    getNEntries() { return fVectorTotal;}
  void     rewindTree()  { fVectorIndex = 0; }
  void     setRewind(Bool_t b) { fRewind = b;}

  void     setVectorIndex(Int_t i) { fVectorIndex = i; }

private:

  TChain  fNuFluxChain;
  TString fKey, fCurrentFileName, fFileName;
  Int_t   fDetector, fNuType;
  Int_t   fVectorIndex, fVectorTotal, fNuFluxEntry;
  Bool_t  fRewind;
  
  Int_t     nChainFluxFiles;
  Int_t     nHaddFluxFiles;

  // ntuple variables (10a)
  Float_t fEnu,  fppi, fxpi[3], fnpi[3], fcospibm, fnorm;
  Float_t fppi0, fxpi0[3], fnpi0[3], fcospi0bm, frnu, fxnu, fynu, fnnu[3];   
  Int_t   fppid, fmode, fnvtx0, fidfd;
  Int_t   fgipart;
  Float_t fgpos0[3], fgvec0[3], fgamom0;

  TBranch *b_Enu, *b_ppid, *b_mode, *b_ppi, *b_xpi; 
  TBranch *b_npi, *b_cospibm, *b_norm;
  TBranch *b_nvtx0, *b_ppi0, *b_xpi0, *b_npi0, *b_cospi0bm; 
  TBranch *b_rnu, *b_xnu, *b_ynu;
  TBranch *b_nnu, *b_idfd;  
  TBranch *b_gipart, *b_gpos0, *b_gvec0, *b_gamom0;

  // 10d
  Int_t   fng;
  Float_t fgpx[12], fgpy[12], fgpz[12], fgcosbm[12];
  Float_t fgvx[12], fgvy[12], fgvz[12];
  Int_t   fgpid[12], fgmec[12]; 
  Float_t fEnusk, fnormsk, fanorm;

  TBranch *b_ng;
  TBranch *b_gpx, *b_gpy, *b_gpz, *b_gcosbm;
  TBranch *b_gvx, *b_gvy, *b_gvz;
  TBranch *b_gpid, *b_gmec; 
  TBranch *b_Enusk, *b_normsk, *b_anorm; 


  // 11a
  Int_t           fgmat[12];
  Float_t         fgdistc[12];    //[ng]       -> NuGdistc
  Float_t         fgdistal[12];   //[ng]       -> NuGdistal
  Float_t         fgdistti[12];   //[ng]       -> NuGdistti
  Float_t         fgdistfe[12];   //[ng]       -> NuGdistfe

  TBranch *b_gmat;
  TBranch *b_gdistc;
  TBranch *b_gdistal;
  TBranch *b_gdistti;
  TBranch *b_gdistfe;
  
  // version tree variables (10d)
  TChain   fVerTree;
  TString  fVerKey;
  TString  fNDKey;
  Int_t    fFileIndex;
  Float_t        fversion, version;
  Int_t          ftuneid, tuneid;
  Int_t          fntrig;
  Int_t          fpint;
  Float_t        fbpos[2];
  Float_t        fbtilt[2];
  Float_t        fbrms[2];
  Float_t        femit[2];
  Float_t        falpha[2];
  Float_t        fhcur[3];
  Int_t          frand;
  Int_t          frseed[2];
   
  TBranch        *b_version;   //!
  TBranch        *b_tuneid;   //!
  TBranch        *b_ntrig;   //!
  TBranch        *b_pint;   //!
  TBranch        *b_bpos;   //!
  TBranch        *b_btilt;   //!
  TBranch        *b_brms;   //!
  TBranch        *b_emit;   //!
  TBranch        *b_alpha;   //!
  TBranch        *b_hcur;   //!
  TBranch        *b_rand;   //!
  TBranch        *b_rseed;   //! 

  // ND information (from h3000 tree)

  // Flux plane position in ND280 coordinates
  Float_t        fFluxPlanePos[99][3];

};
#endif
