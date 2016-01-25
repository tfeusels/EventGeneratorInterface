//////////////////////////////////////////////////////////////////////////
// 
// TNDelements.h
// Class for obtaining elements in nd280 geometry 
//
//////////////////////////////////////////////////////////////////////////

#ifndef TNDelements_h
#define TNDelements_h

#include "TMath.h"
#include "TGeoManager.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "PeriodicTable.h"
#include "TTree.h"

const Int_t kAmax = 300;

class TNDelements {

public:

  // Get elements in full geometry
  TNDelements(TGeoManager *aTheGeometry){

    theGeometry = aTheGeometry;

    hZA = new TH2D("hZA","hZA",kZmax,0,kZmax,kAmax,0,kAmax);

    TIter next(theGeometry->GetListOfMaterials());   
    TGeoMixture *currMixture;
    while ((currMixture=(TGeoMixture*)next())) {

      Double_t *Z = currMixture->GetZmixt();      
      Double_t *A = currMixture->GetAmixt();
      Double_t *W = currMixture->GetWmixt();

      for (int ie=0; ie<currMixture->GetNelements(); ie++) {

	// Get current largest weighted isotope
	Double_t W_temp=0;
	for (int A_bin=1; A_bin<=kAmax; A_bin++) {
	  W_temp = hZA->GetBinContent(TMath::Nint(Z[ie])+1,A_bin);
	  if (W_temp) break;
	}
	
	// Clear Z column if this isotope has higher abundance
	// and store it
	if (W[ie]>W_temp) {
	  for (int A_bin=1; A_bin<=kAmax; A_bin++) {
	    hZA->SetBinContent(TMath::Nint(Z[ie])+1,A_bin,0);
	  }
	  hZA->Fill(TMath::Nint(Z[ie]),TMath::Nint(A[ie]),W[ie]);
	  
	}
	
      }
      
    }


    // Now count how many elements were found
    fNelements=0;
    for (int Z_bin=1; Z_bin<=kZmax; Z_bin++) {
      for (int A_bin=1; A_bin<=kAmax; A_bin++) {
	if (!hZA->GetBinContent(Z_bin,A_bin)) continue;
	
	fZ[fNelements] = Z_bin-1;
	fA[fNelements] = A_bin-1;
	fNelements++;
      }
    }
    
  }



  // Get elements from numc file
  TNDelements(TTree *aTheNumcChain){

    Double_t fgd1_fv[2][3] = {-0.800, -0.745, 0.140,
			      0.800,  0.855 , 0.444 };
    Double_t fgd2_fv[2][3] = {-0.800, -0.745, 1.530,
			      0.800,  0.855 , 1.803 };

    theNumcChain = aTheNumcChain;

    Double_t        EvtVtx[4];
    TBranch        *b_EvtVtx;   //!
   
    Int_t           StdHepPdg[80];   //[StdHepN]
    TBranch        *b_StdHepPdg;   //!

    theNumcChain->SetBranchAddress("EvtVtx", EvtVtx, &b_EvtVtx);
    theNumcChain->SetBranchAddress("StdHepPdg", StdHepPdg, &b_StdHepPdg);

    hZA = new TH2D("hZA","hZA",kZmax,0,kZmax,kAmax,0,kAmax);

    Long64_t nentries = theNumcChain->GetEntriesFast();

    for (Long64_t ientry=0; ientry<nentries;ientry++) {

      if (!(ientry%100000)) std::cout << "Processing event " << ientry << " / " << nentries << std::endl;

      theNumcChain->GetEntry(ientry);
   
      // In FV
      Int_t inFGD1=1;
      Int_t inFGD2=1;
      for (int i=0; i<3; i++) {
	inFGD1 = (fgd1_fv[0][i] < EvtVtx[i] && EvtVtx[i] < fgd1_fv[1][i]) && inFGD1;
	inFGD2 = (fgd2_fv[0][i] < EvtVtx[i] && EvtVtx[i] < fgd2_fv[1][i]) && inFGD2;
      }
      //if (!(inFGD1 || inFGD2)) continue;

      // PDG nuclear code format
      int nucID = abs(StdHepPdg[1]) - 1000000000;
      int nStrange = nucID/10000000;
      int Z = (nucID - nStrange*10000000)/10000;
      int A = (nucID - nStrange*10000000 - Z*10000)/10;
      int I = (nucID - nStrange*10000000 - Z*10000 - A*10);

      hZA->Fill(Z,A);
      
    }
    
    // Now get most abundant isotope for each element
    fNelements=0;
    for (int Z_bin=1; Z_bin<=kZmax; Z_bin++) {

      fW[fNelements] = 0;
      Double_t max_A = 0;

      for (int A_bin=1; A_bin<=kAmax; A_bin++) {
	Double_t temp_A = hZA->GetBinContent(Z_bin,A_bin);
	
	if (!temp_A) continue;
	else if (temp_A > max_A) {
	  max_A = temp_A;
	  fZ[fNelements] = Z_bin-1;
	  fA[fNelements] = A_bin-1;
	}
	fW[fNelements] += temp_A;
      }
      if (max_A) {
	fW[fNelements] = fW[fNelements]/hZA->GetEntries();
	labels[fNelements] = element_symbol[fZ[fNelements]-1];
	fNelements++;      
      }
    }

    

    
  }



  // Get elements from composition file
  TNDelements(TH1D *aTheCompHist){	

    theCompHist = aTheCompHist;
  
    hZA = new TH2D("hZA","hZA",kZmax,0,kZmax,kAmax,0,kAmax);

    fNelements=0;
    for (int Z_bin=1; Z_bin<=theCompHist->GetNbinsX(); Z_bin++) {

      Double_t Z_wt = theCompHist->GetBinContent(Z_bin);

      if (Z_wt) {
	fZ[fNelements] = Z_bin-1;
	fA[fNelements] = a_stable[Z_bin-2];
	fW[fNelements] = Z_wt;

	hZA->Fill(Z_bin-1,a_stable[Z_bin-2]);

	fNelements++;      
      }

    }
    
  }








  // Print
  virtual void Print() {
    std::cout << "Z  A   Name   Weight" << std::endl;
    for (int ie=0; ie<fNelements; ie++) {
      std::cout << fZ[ie] << " " << fA[ie] << " " <<  element_symbol[fZ[ie]-1] << " " << element_name[fZ[ie]-1] << " " << fW[ie] << std::endl;
    }
  }


  TH2D* getZAhisto() {return hZA;}

  Int_t  Nelements() {return fNelements;}
  Int_t* Z() {return fZ;}
  Int_t* A() {return fA;}
  Double_t* W() {return fW;}
  char **getLabels() {return labels;}
  
  TString getElementSymbol(int ie) {
    TString sym = element_symbol[fZ[ie]-1];
    return sym;
  }
  TString getElementName(int ie) {
    TString name = element_name[fZ[ie]-1];
    return name;
  }
  
  ~TNDelements() {
    hZA->Delete();
    
  };



private:


  TGeoManager *theGeometry;
  TTree       *theNumcChain;
  TH1D       *theCompHist;
  Int_t    fNelements;
  TH2D    *hZA;
  Int_t    fZ[kZmax], fA[kAmax];
  Double_t fW[kZmax];
  char     *labels[kZmax];
 
};
#endif
