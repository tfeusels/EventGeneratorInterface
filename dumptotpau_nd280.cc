////////////////////////////////////////////////////////////////////////////////////
//
// Purpose:
//    Creates cross section histograms (total and each mode)
//    for all nuclei present in ND280. Can also create an 
//    effective cross section weighted by atomic abundance in
//    FGD1+2 FV.
// 
// Usage Options: 
// -i Can be two types:
//      1) ND280 geometry file: In this case cross sections are output for all 
//                              nuclei present in the geometry.
//      2) ND280 numc file: Select interactions that occur in the FGD1+2 FV
//                          and determines atomic abundance from truth info
//      3) Composition file: Containing histogram named "hZ_comp" with abundance
//                           of each element 
//         
// -f Flux file for calculating flux_x_xsec and determining atomic abundance from 
//    numc file:
//      e.g. http://www.t2k.org/beam/NuFlux/FluxRelease/10a/fluka_flux/nd5_fine
//    In principle this should correspond to the flux used to generated the numc
//
// -o Output file containing all the histograms
//
//    Note you must set RANFILE to run:  export RANFILE=../neutsmpl/random.tbl
//
////////////////////////////////////////////////////////////////////////////////////

#include <iostream>  
#include <cstdlib>
#include "TFile.h"
#include "TGeoManager.h"
#include "TKey.h"
#include "TPie.h"
#include "TNDelements.h"
#include "neutmodesCPP.h"

int getArgs(int argc, char* argv[]);
TString  inFile            = "";
TString  outputFile        = "nd280_xsecs.root";
TString  fluxFile          = "nd5_fluka_fine.root";

extern "C" {
  void neutxs_(int* nu_in,float* e_in,
		    int* z_in, int* a_in, int* mode_in,
		    float* totxs);
}


int main(int argc, char *argv[]){

 // process the arguments
  int args = getArgs(argc, argv);
  if(args != 0){
    std::cerr << "Usage " << std::endl;
    return 0;
  }

  // open the geometry file
  TFile fnd280in(inFile.Data());
  
  // Get TGeoManager key name
  TIter next(fnd280in.GetListOfKeys());
  TKey *key;
  Int_t foundGeoKey=0;
  Int_t foundNumcKey=0;
  Int_t foundCompKey=0;
  TString keyName, className;
  while ((key=(TKey*)next())) {
    keyName   = key->GetName();
    className = key->GetClassName();

    if (!className.CompareTo("TGeoManager") && keyName.BeginsWith("ND280Geometry")) {
      foundGeoKey=1;
      break;
    }
    else if (!className.CompareTo("TTree") && keyName.BeginsWith("nRooTracker")) {
      foundNumcKey=1;
      break;
    }
    else if (!className.CompareTo("TH1D") && keyName.BeginsWith("hZ_comp")) {
      foundCompKey=1;
      break;
    }
  }
  if (!(foundGeoKey || foundNumcKey || foundCompKey)) {
    std::cerr << "Error: Could not find TGeoManager, TTree or TH1D object in " << inFile << std::endl;
    exit (-1);
  }

  std::cout << "Input file = " << inFile  << std::endl;
  std::cout << "Object class = " << className  << std::endl;
  std::cout << "Key Name = " << keyName  << std::endl;

  // determine elements
  TGeoManager *geom; 
  TTree *tree;
  TNDelements *nd280elements;
  TH1D *hZ_comp;
  if (foundGeoKey) {
    geom = (TGeoManager*)fnd280in.Get(keyName);
    nd280elements = new TNDelements(geom);
  }
  else if (foundNumcKey) {
    tree = (TTree*)fnd280in.Get(keyName);
    nd280elements = new TNDelements(tree);
  }
  else if (foundCompKey) {
    hZ_comp = new TH1D(*(TH1D*)fnd280in.Get(keyName));
    nd280elements = new TNDelements(hZ_comp);
  }
  
  nd280elements->Print();


  const int kNelements = nd280elements->Nelements();
  Int_t *Z = nd280elements->Z();
  Int_t *A = nd280elements->A();

  // Now make cross section histograms
  Double_t POT = 0.0294;  // x10^21 POT
  Int_t  enuLow = 0;  // GeV
  Int_t enuHigh = 10;
  Double_t binWidth = 0.05;
  Int_t nBins = (enuHigh-enuLow)/binWidth;

  char  nuName[4][6]  = {"numu","nue","numub","nueb"};
  char nuTitle[4][13] = {"#nu#mu","#nue","#bar{#nu#mu}","#bar{#nue}"};
  int     nuID[4]    = {14,12,-14,-12};
 

  TFile *theOutputFile = new TFile(outputFile,"RECREATE");
  
  // Store array of elements as histogram
  TH1I *hZ = new TH1I("hZ","hZ;Index (ie);Z",kNelements,0,kNelements);
  for (int ie=0; ie<kNelements; ie++)
    hZ->SetBinContent(ie+1,Z[ie]);

  (nd280elements->getZAhisto())->Write();

  TH1D *h_xsec[kNelements][4][kNneutModes];

  for (int ie=0; ie<kNelements; ie++) {

    std::cout << "Generating Xsecs for " << element_name[Z[ie]-1] << ", Z=" << Z[ie] << ", A=" << A[ie] << std::endl;

    for (int inu=0; inu<4; inu++) {
      
      for (int imode=0; imode<kNneutModes; imode++) {

	h_xsec[ie][inu][imode] = new TH1D(Form("%s_xsec_%s_%s",element_symbol[Z[ie]-1],nuName[inu],neutModeName[imode]),Form("^{%d}_{%d}%s %s;True E_{%s} [GeV];#sigma [cm^{2}/atom]",A[ie],Z[ie],element_symbol[Z[ie]-1],neutModeTitle[imode],nuTitle[inu]),nBins,enuLow,enuHigh);
	
	for (float enu=enuLow+binWidth/2; enu<enuHigh; enu+=binWidth) {
	  float neutxs;
	  neutxs_(&nuID[inu], &enu, &Z[ie], &A[ie], &neutModeID[imode], &neutxs);
	  h_xsec[ie][inu][imode]->Fill(enu,neutxs*pow(10,-38));
	} 
      }
    }
  }
  std::cout << std::endl;

  
  // Following section is for calculating effective cross section 
  TFile *theFluxFile = new TFile(fluxFile);
  TH1D *h_flux[4];
  TH1D *h_flux_x_xsec[kNelements][5];

  theOutputFile->cd();
  for (int inu=0; inu<4; inu++) {
    h_flux[inu] = new TH1D(*(TH1D*)theFluxFile->Get(Form("enu_nd5_10a_norm_%s",nuName[inu])));

    h_flux[inu]->Scale(POT*h_flux[inu]->GetBinWidth(1));
    h_flux[inu]->GetYaxis()->SetTitle("Flux [/cm^{2}]");


    for (int ie=0; ie<kNelements; ie++) {
      h_flux_x_xsec[ie][inu] = new TH1D(*h_flux[inu]);
      h_flux_x_xsec[ie][inu]->SetName(Form("%s_x_xsec_%s",nuName[inu],element_symbol[Z[ie]-1]));
      h_flux_x_xsec[ie][inu]->SetTitle(Form("%s Flux #times #sigma_{%s}",nuName[inu],element_symbol[Z[ie]-1]));
      h_flux_x_xsec[ie][inu]->GetYaxis()->SetTitle(Form("%s Flux #times #sigma_{%s}",nuTitle[inu],element_symbol[Z[ie]-1]));
      h_flux_x_xsec[ie][inu]->Multiply(h_xsec[ie][inu][0]);
    }

  }

  // Sum flux_x_xsec
  for (int ie=0; ie<kNelements; ie++) {
    h_flux_x_xsec[ie][4] = new TH1D(*h_flux_x_xsec[ie][0]);
    h_flux_x_xsec[ie][4]->Reset();
    h_flux_x_xsec[ie][4]->SetName(Form("sum_flux_x_xsec_%s",element_symbol[Z[ie]-1]));
    h_flux_x_xsec[ie][4]->GetYaxis()->SetTitle(Form("Flux #times #sigma_{%s}",element_symbol[Z[ie]-1]));
    h_flux_x_xsec[ie][4]->SetTitle(Form("Sum of Flux #times #sigma_{%s}",element_symbol[Z[ie]-1]));

    for (int inu=0; inu<4; inu++)
      h_flux_x_xsec[ie][4]->Add(h_flux_x_xsec[ie][inu]);
  }	
  /**/

  
  // Now do some calculations to find molecule ratio
  Double_t *real_ratio;

  // Pie Chart colours
  Int_t pieColors[kNelements];
  for (int ie=0; ie<kNelements; ie++) 
    pieColors[ie] = ie+2;

  char *pieLabels[kNelements];
    

  if (foundNumcKey) {

    Double_t *W = nd280elements->W();
   

    // Pie Chart labels
    for (int ie=0; ie<kNelements; ie++) 
      pieLabels[ie] = Form("%s (%.02g\%)",(nd280elements->getLabels())[ie],W[ie]*100);
    

    TPie *pZN = new TPie("pZN","pZN",kNelements,W,pieColors,pieLabels);
    pZN->Write();
    

    Double_t ratio[kNelements][kNelements];      
    Double_t numEvents[kNelements];
    for (int ie=0; ie<kNelements; ie++) {
      numEvents[ie] = h_flux_x_xsec[ie][4]->Integral();
    }
    for (int ie_num=0; ie_num<kNelements; ie_num++) {
      for (int ie_den=0; ie_den<kNelements; ie_den++) {

	ratio[ie_num][ie_den] = (W[ie_num]/W[ie_den])*(numEvents[ie_den]/numEvents[ie_num]);
      }
    }

    real_ratio = new Double_t[kNelements];
    
    std::cout << "Element Nint(%) Atom Abundance" << std::endl;
    
    for (int ie_num=0; ie_num<kNelements; ie_num++) {
      for (int ie_den=0; ie_den<kNelements; ie_den++) {
	real_ratio[ie_num] += 1./ratio[ie_num][ie_den];
      }
      real_ratio[ie_num] = 1./real_ratio[ie_num];
    
      std::cout << element_name[Z[ie_num]-1] << " " << W[ie_num] << " " << real_ratio[ie_num] << std::endl;
    }




  }
  else if (foundCompKey) {
    real_ratio = nd280elements->W();
  }

   
  if (foundCompKey || foundNumcKey) {

    for (int ie=0; ie<kNelements; ie++) 
      pieLabels[ie] = Form("%s (%.02g\%)",(nd280elements->getLabels())[ie],real_ratio[ie]*100);
    TPie *pZW = new TPie("pZW","pZW",kNelements,real_ratio,pieColors,pieLabels);
    pZW->Write();


    // Now make effective cross sections
    TH1D *eff_xsec[4][kNneutModes];
    for (int inu=0; inu<4; inu++) {
      for (int imode=0; imode<kNneutModes; imode++) {
      
	eff_xsec[inu][imode] = new TH1D(Form("eff_xsec_%s_%s",nuName[inu],neutModeName[imode]),Form("FGD1+2 FV Eff. #sigma (%s);True E_{%s} [GeV];#sigma [cm^{2}/atom]",neutModeTitle[imode],nuTitle[inu]),nBins,enuLow,enuHigh);

	for (int ie=0; ie<kNelements; ie++) {
	  h_xsec[ie][inu][imode]->Scale(real_ratio[ie]);
	  eff_xsec[inu][imode]->Add(h_xsec[ie][inu][imode]);
	  h_xsec[ie][inu][imode]->Scale(1./(real_ratio[ie]));  // Back to normal for output
	}
      }
    }
  }
  std::cout << std::endl;


  theOutputFile->Write();

  std::cout << "Output file = " << outputFile.Data() << std::endl;

  theOutputFile->Close();

  return 0;
  
}




int getArgs(int argc, char* argv[]){

  while( (argc > 1) && (argv[1][0] == '-') ){

    switch(argv[1][1]){

    case 'i': 
      inFile = argv[2];
      ++argv; --argc;
      break;


    case 'o': 
      outputFile = argv[2];
      ++argv; --argc;
      break;

    case 'f': 
      fluxFile = argv[2];
      ++argv; --argc;
      break;
    }
    ++argv; --argc;
  }
  return 0;
  
}

