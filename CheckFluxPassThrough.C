//
// This macro is for checking internal consistency of flux pass through information
// and also to check GENIE and NEUT are passing the same info and in the same format. 
// 
// Loops over a rootracker (gRooTracker or nRooTracker) file and compares pass-through
// info to the original data in the input flux file used for vector generation. This
// needs to be run from the directory containing the original flus file:
// [root].x CheckFluxPassThrough.C("path/to/rootracker/file.root") 
//
// J.Dobson - 4th Feb 2011
//
#include <iostream>

using std::cout;
using std::endl;

// A couple of functions for comparing arrays of int's and floats
double fASmallNum = 0.000005;
bool Compare(float a, float b); 
bool Compare(int a, int b);
bool Compare(float a[], float b[], int ni);
bool Compare(int a[], int b[], int ni);
// Convert back to Geant code to allow comparison to flux id's
int PdgToGeant(int pdg);

//_____________________________________________________________________________
//
// Main function for checking that have passed-through all the JNuBeam information
// correctly. Loops over all rootracker events and compares values of pass-through
// information to the corresponding original flux information.
void CheckFluxPassThrough(char * rootracker_filename){

  // Look for a gRooTracker or nRooTracker input file
  TFile * infile = new TFile(rootracker_filename); 
  assert(infile);
  TTree * rooTrackerTree = (TTree*)infile->Get("gRooTracker");
  if(!rooTrackerTree){
    rooTrackerTree = (TTree*)infile->Get("nRooTracker");
    if(!rooTrackerTree){
      printf("Could not find gRooTracker or nRooTracker trees!\n"); 
      exit(1);
    }
  }
  assert(rooTrackerTree->GetEntries()>0);

  // Get the original flux file name
  TObjString* NuFileName; rooTrackerTree->SetBranchAddress("NuFileName", &NuFileName);       
  assert(rooTrackerTree->GetEntry(0));
  TString fluxfilename = NuFileName->GetString();
  fluxfilename = fluxfilename(0, fluxfilename.Index(":"));

  printf("\nFound tree: %s \n", rooTrackerTree->GetName());;
  printf("Checking against flux file name: %s \n\n", fluxfilename.Data());
  // Load the two input trees expected in 11a flux
  TChain * h1000_ch = new TChain("h1000");
  h1000_ch->Add(fluxfilename);
  assert(h1000_ch->GetEntries() > 0);
  TChain * h3002_ch = new TChain("h3002");
  h3002_ch->Add(fluxfilename);
  assert(h3002_ch->GetEntries() > 0);
 
  // Generate some classes for the original flux trees
  h3002_ch->MakeClass("MCh3002"); 
  h1000_ch->MakeClass("MCh1000"); 
  // Load the make tree projects and make some instances
  gROOT->ProcessLine(".L MCh1000.C"); 
  gROOT->ProcessLine(".L MCh3002.C");  
  MCh3002  h3002;  
  MCh1000  h1000;  
 
  // Set branch addresses of input gRooTracker file to compare to the original
  // flux file variables. Do this manually as want to catch any differences in the
  // way NEUT and GENIE are saving them, i.e. any typos would show up this way.
  Long64_t    NuFluxEntry; rooTrackerTree->SetBranchAddress("NuFluxEntry", &NuFluxEntry);       
  int         NuParentPdg; rooTrackerTree->SetBranchAddress("NuParentPdg", &NuParentPdg);     
  int         NuParentDecMode; rooTrackerTree->SetBranchAddress("NuParentDecMode", &NuParentDecMode);
  double       NuParentDecP4 [4]; rooTrackerTree->SetBranchAddress("NuParentDecP4", NuParentDecP4);
  double       NuParentDecX4 [4]; rooTrackerTree->SetBranchAddress("NuParentDecX4", NuParentDecX4);
  float       NuCospibm; rooTrackerTree->SetBranchAddress("NuCospibm", &NuCospibm);
  float       NuNorm; rooTrackerTree->SetBranchAddress("NuNorm", &NuNorm);
  double       NuParentProP4 [4]; rooTrackerTree->SetBranchAddress("NuParentProP4", NuParentProP4);
  double       NuParentProX4 [4]; rooTrackerTree->SetBranchAddress("NuParentProX4", NuParentProX4);
  float       NuCospi0bm; rooTrackerTree->SetBranchAddress("NuCospi0bm", &NuCospi0bm);
  float       NuRnu; rooTrackerTree->SetBranchAddress("NuRnu", &NuRnu);
  float       NuXnu [2]; rooTrackerTree->SetBranchAddress("NuXnu", NuXnu);
  int         NuIdfd; rooTrackerTree->SetBranchAddress("NuIdfd", &NuIdfd);
  int         NuGipart; rooTrackerTree->SetBranchAddress("NuGipart", &NuGipart);
  float       NuGpos0[3]; rooTrackerTree->SetBranchAddress("NuGpos0", NuGpos0);
  float       NuGvec0[3]; rooTrackerTree->SetBranchAddress("NuGvec0", NuGvec0);
  float       NuGamom0; rooTrackerTree->SetBranchAddress("NuGamom0", &NuGamom0); 
  int         NuNg; rooTrackerTree->SetBranchAddress("NuNg", &NuNg);     
  float       NuGp[12][3]; rooTrackerTree->SetBranchAddress("NuGp", NuGp);
  float       NuGv[12][3]; rooTrackerTree->SetBranchAddress("NuGv", NuGv);
  float       NuGcosbm[12]; rooTrackerTree->SetBranchAddress("NuGcosbm", NuGcosbm);
  int         NuGmec[12]; rooTrackerTree->SetBranchAddress("NuGmec", NuGmec);
  float       NuEnusk; rooTrackerTree->SetBranchAddress("NuEnusk", &NuEnusk);
  float       NuNormsk; rooTrackerTree->SetBranchAddress("NuNormsk", &NuNormsk);
  float       NuAnorm; rooTrackerTree->SetBranchAddress("NuAnorm", &NuAnorm);
  int         NuGmat[12]; rooTrackerTree->SetBranchAddress("NuGmat", NuGmat);
  float       NuGdistc[12]; rooTrackerTree->SetBranchAddress("NuGdistc", NuGdistc);
  float       NuGdistal[12]; rooTrackerTree->SetBranchAddress("NuGdistal", NuGdistal);
  float       NuGdistti[12]; rooTrackerTree->SetBranchAddress("NuGdistti", NuGdistti);
  float       NuGdistfe[12]; rooTrackerTree->SetBranchAddress("NuGdistfe", NuGdistfe);
  float       NuVersion; rooTrackerTree->SetBranchAddress("NuVersion", &NuVersion);
  int         NuTuneid; rooTrackerTree->SetBranchAddress("NuTuneid", &NuTuneid);
  int         NuNtrig; rooTrackerTree->SetBranchAddress("NuNtrig", &NuNtrig);
  int         NuPint; rooTrackerTree->SetBranchAddress("NuPint", &NuPint);
  float       NuBpos[2]; rooTrackerTree->SetBranchAddress("NuBpos", NuBpos);
  float       NuBtilt[2]; rooTrackerTree->SetBranchAddress("NuBtilt", NuBtilt);
  float       NuBrms[2]; rooTrackerTree->SetBranchAddress("NuBrms", NuBrms);
  float       NuEmit[2]; rooTrackerTree->SetBranchAddress("NuEmit", NuEmit);
  float       NuAlpha[2]; rooTrackerTree->SetBranchAddress("NuAlpha", NuAlpha);
  float       NuHcur[3]; rooTrackerTree->SetBranchAddress("NuHcur", NuHcur);
  int         NuRand; rooTrackerTree->SetBranchAddress("NuRand", &NuRand);
  int         NuRseed[2]; rooTrackerTree->SetBranchAddress("NuRseed", NuRseed);

  // Now loop over all entries in the rootracker file and check that the values stored for
  // for the flux pass-through are the same as those in the original tree.
  int global_counter = 0;
  printf("\nLooping over all rootracker entries and checking pass-through info:\n");
  for(int i = 0; i< rooTrackerTree->GetEntries(); i++){
    rooTrackerTree->GetEntry(i);
    h3002.GetEntry(NuFluxEntry);
    h1000.GetEntry(0);

    if(i%500 ==0){printf("   -> checked %g\%\n", 100.0*i/rooTrackerTree->GetEntries());}
    
    int nfailed = 0;
    if( !Compare(NuParentDecP4[0], (h3002.npi[0]*h3002.ppi)) ){ nfailed++; printf("Failed check for: NuParentDecP4[0]!\n"); }
    if( !Compare(NuParentDecP4[1], (h3002.npi[1]*h3002.ppi)) ){ nfailed++; printf("Failed check for: NuParentDecP4[1]!\n"); }
    if( !Compare(NuParentDecP4[2], (h3002.npi[2]*h3002.ppi)) ){ nfailed++; printf("Failed check for: NuParentDecP4[2]!\n"); }
    if( !Compare(NuParentDecX4[0], (h3002.xpi[0])) ){ nfailed++; printf("Failed check for: NuParentDecX4[0]!\n"); }
    if( !Compare(NuParentDecX4[1], (h3002.xpi[1])) ){ nfailed++; printf("Failed check for: NuParentDecX4[1]!\n"); }
    if( !Compare(NuParentDecX4[2], (h3002.xpi[2])) ){ nfailed++; printf("Failed check for: NuParentDecX4[2]!\n"); }
    for(int j=0;j<NuNg;j++){
      if( !Compare(NuGv[j][0], h3002.gvx[j]) ){ nfailed++; printf("Failed check for: NuGv[j][0]!\n"); }
      if( !Compare(NuGv[j][1], h3002.gvy[j]) ){ nfailed++; printf("Failed check for: NuGv[j][1]!\n"); }
      if( !Compare(NuGv[j][2], h3002.gvz[j]) ){ nfailed++; printf("Failed check for: NuGv[j][2]!\n"); }
      if( !Compare(NuGp[j][0], h3002.gpx[j]) ){ nfailed++; printf("Failed check for: NuGp[j][0]!\n"); }
      if( !Compare(NuGp[j][1], h3002.gpy[j]) ){ nfailed++; printf("Failed check for: NuGp[j][1]!\n"); }
      if( !Compare(NuGp[j][2], h3002.gpz[j]) ){ nfailed++; printf("Failed check for: NuGp[j][2]!\n"); }
    } 
    if( !Compare( PdgToGeant(NuParentPdg), h3002.ppid) ){ nfailed++; printf("Failed check for: NuParentPdg!\n");}
    if( !Compare(NuParentDecMode, h3002.mode) ){ nfailed++; printf("Failed check for: NuParentDecMode!\n"); }
    if( !Compare(NuCospibm, h3002.cospibm) ){ nfailed++; printf("Failed check for: NuCospibm!\n"); }
    if( !Compare(NuNorm, h3002.norm) ){ nfailed++; printf("Failed check for: NuNorm!\n"); } 
    if( !Compare(NuCospi0bm, h3002.cospi0bm) ){ nfailed++; printf("Failed check for: NuCospi0bm!\n"); }
    if( !Compare(NuRnu, h3002.rnu) ){ nfailed++; printf("Failed check for: NuRnu!\n"); }
    if( !Compare(NuXnu[0], h3002.xnu) ){ nfailed++; printf("Failed check for: NuXnu[0]!\n"); }
    if( !Compare(NuXnu[1], h3002.ynu) ){ nfailed++; printf("Failed check for: NuXnu[1]!\n");}
    if( !Compare(NuIdfd, h3002.idfd) ){ nfailed++; printf("Failed check for: NuIdfd!\n"); }
    if( !Compare(NuGipart, h3002.gipart) ){ nfailed++; printf("Failed check for: NuGipart!\n");}
    if( !Compare(NuGpos0, h3002.gpos0, 3) ){ nfailed++; printf("Failed check for: NuGpos0!\n"); }
    if( !Compare(NuGvec0, h3002.gvec0, 3) ){ nfailed++; printf("Failed check for: NuGvec0!\n"); }
    if( !Compare(NuGamom0, h3002.gamom0) ){ nfailed++; printf("Failed check for: NuGamom0!\n"); }
    if( !Compare(NuNg, h3002.ng) ){ nfailed++; printf("Failed check for: NuNg!\n"); }
    if( !Compare(NuGcosbm, h3002.gcosbm, NuNg) ){ nfailed++; printf("Failed check for: NuGcosbm!\n");}
    if( !Compare(NuGmec, h3002.gmec, NuNg) ){ nfailed++; printf("Failed check for: NuGmec!\n"); }
    if( !Compare(NuEnusk, h3002.Enusk) ){ nfailed++; printf("Failed check for: NuEnusk!\n"); }
    if( !Compare(NuNormsk, h3002.normsk) ){ nfailed++; printf("Failed check for: NuNormsk\n");}
    if( !Compare(NuAnorm, h3002.anorm) ){ nfailed++; printf("Failed check for: NuAnorm!\n"); }
    if( !Compare(NuGmat, h3002.gmat, NuNg) ){ nfailed++; printf("Failed check for: NuGmat!\n"); }
    if( !Compare(NuGdistc, h3002.gdistc, NuNg) ){ nfailed++; printf("Failed check for: NuGdistc!\n");}
    if( !Compare(NuGdistal, h3002.gdistal, NuNg) ){ nfailed++; printf("Failed check for: NuGdistal!\n"); }
    if( !Compare(NuGdistti, h3002.gdistti, NuNg) ){ nfailed++; printf("Failed check for: NuGdistti!\n"); }
    if( !Compare(NuGdistfe, h3002.gdistfe, NuNg) ){ nfailed++; printf("Failed check for: NuGdistfe!\n"); }
    if( !Compare(NuVersion, h1000.version) ){ nfailed++; printf("Failed check for: NuVersion!\n"); }
    if( !Compare(NuTuneid, h1000.tuneid) ){ nfailed++; printf("Failed check for: NuTuneid!\n"); }
    if( !Compare(NuNtrig, h1000.ntrig) ){ nfailed++; printf("Failed check for: NuNtrig!\n"); }
    if( !Compare(NuPint, h1000.pint) ){ nfailed++; printf("Failed check for: NuPint!\n"); }
    if( !Compare(NuBpos, h1000.bpos, 2) ){ nfailed++; printf("Failed check for: NuBpos!\n");}
    if( !Compare(NuBtilt, h1000.btilt, 2) ){ nfailed++; printf("Failed check for: NuBtilt!\n"); }
    if( !Compare(NuBrms, h1000.brms, 2) ){ nfailed++; printf("Failed check for: NuBrms!\n"); }
    if( !Compare(NuEmit, h1000.emit, 2) ){ nfailed++; printf("Failed check for: NuEmit!\n"); }
    if( !Compare(NuAlpha, h1000.alpha, 2) ){ nfailed++; printf("Failed check for: NuAlpha!\n"); }
    if( !Compare(NuHcur, h1000.hcur, 3) ){ nfailed++; printf("Failed check for: NuHcur!\n"); }
    if( !Compare(NuRand, h1000.rand) ){ nfailed++; printf("Failed check for: NuRand!\n"); }
    //if( !Compare(NuRseed, h1000.rseed, 2) ){ nfailed++; printf("Failed check for: NuRseed!\n");}
    if(nfailed > 0){ 
      printf("Have had %d failures for entry %d!\n", nfailed , i ); 
      global_counter++;
    }   
  } // loop over rootracker events
  printf("There were %d events with incorrect pass-through info\n", global_counter);
  if(global_counter == 0){
    printf("\n\n------------ PASSED TEST ------------\n\n");
    return;
  } 
  printf("\n\n------------ PASSED TEST ------------\n\n");
  return;
}

//_____________________________________________________________________________
int PdgToGeant(int pdg){
  switch (pdg){
    case (11): return   3;// return kPdgElectron;     //    11 / e-
    case (-11): return   2;// return kPdgPositron;     //   -11 / e+
    case (13): return   6;// return kPdgMuon;         //    13 / mu-
    case (-13): return   5;// return kPdgAntiMuon;     //   -13 / mu+             
    case (15): return  34;// return kPdgTau;          //    15 / tau-
    case (-15): return  33;// return kPdgAntiTau;      //   -15 / tau+              
    case (211): return   8;// return kPdgPiP;          //   211 / pi+
    case (-211): return   9;// return kPdgPiM;          //  -211 / pi-
    case (111): return   7;// return kPdgPi0;          //   111 / pi0
    case (221): return  17;// return kPdgEta;          //   221 / eta
    case (321): return  11;// return kPdgKP;           //   321 / K+
    case (-321): return  12;// return kPdgKM;           //  -321 / K-
    case (130): return  10;// return kPdgK0L;          //   130 / K0_{long}
    case (310): return  16;// return kPdgK0S;          //   310 / K0_{short}
    case (411): return  35;// return kPdgDP;           //   411 / D+
    case (-411): return  36;// return kPdgDM;           //  -411 / D-
    case (421): return  37;// return kPdgD0;           //   421 / D0
    case (-421): return  38;// return kPdgAntiD0;       //  -421 / \bar{D0}
    case (431): return  39;// return kPdgDPs;          //   431 / D+_{s}
    case (-431): return  40;// return kPdgDMs;          //  -431 / D-_{s}
    case (22): return   1;// return kPdgGamma;        //    22 / photon
    case (23): return  44;// return kPdgZ0;           //    23 / Z
    case (24): return  42;// return kPdgWP;           //    24 / W+
    case (-24): return  43;// return kPdgWM;           //   -24 / W-
    case (2212): return  14;// return kPdgProton;       //  2212
    case (-2212): return  15;// return kPdgAntiProton;   // -2212
    case (2112): return  13;// return kPdgNeutron;      //  2112
    case (-2112): return  25;// return kPdgAntiNeutron;  // -2112
    case (3122): return  18;// return kPdgLambda;       //  3122 / Lambda
    case (-3122): return  26;// return kPdgAntiLambda;   // -3122 / \bar{Lambda}
    case (3222): return  19;// return kPdgSigmaP;       //  3222 / Sigma+
    case (3212): return  20;// return kPdgSigma0;       //  3212 / Sigma0
    case (3112): return  21;// return kPdgSigmaM;       //  3112 / Sigma-
    case (-3112): return  29;// return kPdgAntiSigmaP;   // -3112 / \bar{Sigma+}
    case (-3212): return  28;// return kPdgAntiSigma0;   // -3212 / \bar{Sigma0}
    case (-3112): return  27;// return kPdgAntiSigmaM;   // -3112 / \barSigma-}
    case (3322): return  22;// return kPdgXi0;          //  3322 / Xi0
    case (3312): return  23;// return kPdgXiM;          //  3312 / Xi-
    case (-3322): return  30;// return kPdgAntiXi0;      // -3322 / \bar{Xi0}
    case (-3312): return  31;// return kPdgAntiXiP;      // -3312 / \bar{Xi+}
    case (3334): return  24;// return kPdgOmegaM;       //  3334 / Omega-
    case (-3334): return  32;// return kPdgAntiOmegaP;   // -3334 / \bar{Omega+}
    default: return 0;
  }
}

//_____________________________________________________________________________
bool Compare(double a, double b){ 
  bool arethesame = TMath::Abs(a-b) < fASmallNum;
  if(!arethesame){ printf("Compared: %g to %g --> diff = %g \n", a ,b , a-b); }
  return arethesame;
}

//_____________________________________________________________________________
bool Compare(float a, float b){ 
  bool arethesame = TMath::Abs(a-b) < fASmallNum;
  if(!arethesame){ printf("Compared: %g to %g --> diff = %g \n", a ,b , a-b );}
  return arethesame;
}

//_____________________________________________________________________________
bool Compare(int a, int b){ 
  bool arethesame = TMath::Abs(a-b) < fASmallNum;
  if(!arethesame){ printf("Compared: %g to %g --> diff = %g \n", a , b , a-b ); }
  return arethesame;
}

//_____________________________________________________________________________
bool Compare(float a[], float b[], int ni){
  for(int i=0; i<ni; i++){ if(Compare(a[i],b[i]) == false){ return false;} }
  return true;
}

//_____________________________________________________________________________
bool Compare(int a[], int b[], int ni){ 
  for(int i=0; i<ni; i++){ if(Compare(a[i],b[i]) == false){ return false;} }
  return true;
}


