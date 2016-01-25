#include "TNeutOutput.h"
#include "TVector3.h"

TNeutOutput::TNeutOutput(TString theFilename, Int_t format, Float_t fluxVersion){
  OpenFile(theFilename);
  SetOutputFormat(format);
  InitTree(fluxVersion);
}
void TNeutOutput::OpenFile(TString theFilename){
  fOutputFile = new TFile(theFilename.Data(), "RECREATE");
}

void TNeutOutput::WriteTree(){
  fOutputFile->cd();
  fOutputTree->Write();
}

void TNeutOutput::CloseFile(){
  fOutputFile->cd();
  fOutputFile->Write();
  fOutputFile->Close();
}

void TNeutOutput::InitTree(Float_t fluxVersion){

  if(fFormat == 0){
    fOutputTree = new TTree("h10","neut");
    fOutputTree->Branch("Nev"      ,&fNev     , "Nev/I"          );
    fOutputTree->Branch("Pos"      ,fPos      , "Pos[3]/F"       );
     
    fOutputTree->Branch("Mode"     ,&fMode    , "Mode/I"         );
    fOutputTree->Branch("Numnu"    ,&fNumnu   , "Numnu/I"        );
    fOutputTree->Branch("Ipnu"     ,fIpnu     , "Ipnu[Numnu]/I"  );
    fOutputTree->Branch("Abspnu"   ,fAbspnu   , "Abspnu[Numnu]/F");
    fOutputTree->Branch("Pnu"      ,fPnu      , "Pnu[Numnu][3]/F");
    fOutputTree->Branch("Npar"     ,&fNpar    , "Npar/I"         );
    fOutputTree->Branch("Ipv"      ,fIpv      , "Ipv[Npar]/I"    );
    fOutputTree->Branch("Icrnv"    ,fIcrnv    , "Icrnv[Npar]/I"  );
    fOutputTree->Branch("Pmomv"    ,fPmomv    , "Pmomv[Npar][3]" );
    
    fOutputTree->Branch("xpi"      ,fxpi      , "xpi[3]/F"  );
    fOutputTree->Branch("npi"      ,fnpi      , "npi[3]/F"  );
    fOutputTree->Branch("cospibm"  ,&fcospibm , "cospibm/F" );
    fOutputTree->Branch("ppi"      ,&fppi     , "ppi/F"     );
    fOutputTree->Branch("ppid"     ,&fppid    , "ppid/I"    );
    fOutputTree->Branch("xpi0"     ,fxpi0     , "xpi0[3]/F" );
    fOutputTree->Branch("npi0"     ,fnpi0     , "npi0[3]/F" );
    fOutputTree->Branch("cospi0bm" ,&fcospi0bm, "cospi0bm/F");
    fOutputTree->Branch("ppi0"     ,&fppi0    , "ppi0/F"    );
  } else if(fFormat == 1){
    std::cout << "in rootracker initialization" << std::endl;

    //"n" for neut
    fOutputTree = new TTree("nRooTracker","RooTracker");
    fEvtFlags = NULL;
    fOutputTree->Branch("EvtFlags", "TBits",     &fEvtFlags, 32000, 1);
    fEvtCode = NULL;
    fOutputTree->Branch("EvtCode",  "TObjString", &fEvtCode, 32000, 1);
    fOutputTree->Branch("EvtNum",         &fEvtNum     ,     "EvtNum/I"                );
    //fOutputTree->Branch("NEibound",       &fNEibound    ,    "NEibound/I"              );
    fOutputTree->Branch("EvtXSec",        &fEvtXSec    ,     "EvtXSec/D"               );
    fOutputTree->Branch("EvtDXSec",       &fEvtDXSec   ,     "EvtDXSec/D"              );
    fOutputTree->Branch("EvtWght",        &fEvtWght    ,     "EvtWght/D"               );
    fOutputTree->Branch("EvtProb",        &fEvtProb    ,     "EvtProb/D"               );
    fOutputTree->Branch("EvtVtx",          fEvtVtx     ,     "EvtVtx[4]/D"             );
    fOutputTree->Branch("StdHepN",        &fStdHepN    ,     "StdHepN/I"               );
    fOutputTree->Branch("StdHepPdg",       fStdHepPdg  ,     "StdHepPdg[StdHepN]/I"    );
    fOutputTree->Branch("StdHepStatus",    fStdHepStatus,    "StdHepStatus[StdHepN]/I" );
    fOutputTree->Branch("StdHepX4",        fStdHepX4,        "StdHepX4[StdHepN][4]/D"  );
    fOutputTree->Branch("StdHepP4",        fStdHepP4,        "StdHepP4[StdHepN][4]/D"  );
    fOutputTree->Branch("StdHepPolz",      fStdHepPolz,      "StdHepPolz[StdHepN][3]/D");
    fOutputTree->Branch("StdHepFd",        fStdHepFd,        "StdHepFd[StdHepN]/I"     );
    fOutputTree->Branch("StdHepLd",        fStdHepLd,        "StdHepLd[StdHepN]/I"     );
    fOutputTree->Branch("StdHepFm",        fStdHepFm,        "StdHepFm[StdHepN]/I"     );
    fOutputTree->Branch("StdHepLm",        fStdHepLm,        "StdHepLm[StdHepN]/I"     );

    // NEUT Native VCWORK block
    fOutputTree->Branch("NEneutmode",     &fNEneutmode    ,  "NEneutmode/I"            );
    fOutputTree->Branch("NEnvc",          &fNEnvc    ,       "NEnvc/I"                 );
    fOutputTree->Branch("NEipvc",          fNEipvc    ,      "NEipvc[NEnvc]/I"         );
    fOutputTree->Branch("NEpvc",           fNEpvc    ,       "NEpvc[NEnvc][3]/F"       );
    fOutputTree->Branch("NEiorgvc",        fNEiorgvc    ,    "NEiorgvc[NEnvc]/I"       );
    fOutputTree->Branch("NEiflgvc",        fNEiflgvc    ,    "NEiflgvc[NEnvc]/I"       ); 
    fOutputTree->Branch("NEicrnvc",        fNEicrnvc    ,    "NEicrnvc[NEnvc]/I"       );
    fOutputTree->Branch("NEcrsx",         &fNEcrsx      ,    "NEcrsx/F"                );
    fOutputTree->Branch("NEcrsy",         &fNEcrsy      ,    "NEcrsy/F"                );
    fOutputTree->Branch("NEcrsz",         &fNEcrsz      ,    "NEcrsz/F"                );
    fOutputTree->Branch("NEcrsphi",       &fNEcrsphi     ,   "NEcrsphi/F"              );


    // FSIHIST block : NEUT version > 5.0.7
    fOutputTree->Branch("NEnvert",        &fNEnvert    ,     "NEnvert/I"               );
    fOutputTree->Branch("NEposvert",       fNEposvert    ,   "NEposvert[NEnvert][3]/F" );
    fOutputTree->Branch("NEiflgvert",      fNEiflgvert    ,  "NEiflgvert[NEnvert]/I"   );
    fOutputTree->Branch("NEnvcvert",      &fNEnvcvert    ,   "NEnvcvert/I"             );
    fOutputTree->Branch("NEdirvert",       fNEdirvert    ,   "NEdirvert[NEnvcvert][3]/F");
    fOutputTree->Branch("NEabspvert",      fNEabspvert    ,  "NEabspvert[NEnvcvert]/F" );
    fOutputTree->Branch("NEabstpvert",     fNEabstpvert   ,  "NEabstpvert[NEnvcvert]/F");
    fOutputTree->Branch("NEipvert",        fNEipvert    ,    "NEipvert[NEnvcvert]/I"   );
    fOutputTree->Branch("NEiverti",        fNEiverti    ,    "NEiverti[NEnvcvert]/I"   );
    fOutputTree->Branch("NEivertf",        fNEivertf    ,    "NEivertf[NEnvcvert]/I"   );

    if(fNuFileName) { 
      std::cout << "Warning in TNeutOutput.cc: fNuFileName has non-zero memory address at initialization" << std::endl;
      fNuFileName = NULL;
    }
    
    fOutputTree->Branch("NuFileName",    "TObjString", &fNuFileName,       1000*32000, 1   );
        
    fOutputTree->Branch("NuFluxEntry",   &fNuFluxEntry,     "NuFluxEntry/L"           );
   
    fOutputTree->Branch("NuParentPdg",    &fNuParentPdg,     "NuParentPdg/I"           );
    fOutputTree->Branch("NuParentDecMode",&fNuParentDecMode, "NuParentDecMode/I"       );
    fOutputTree->Branch("NuParentDecP4",   fNuParentDecP4,   "NuParentDecP4[4]/D"      );
    fOutputTree->Branch("NuParentDecX4",   fNuParentDecX4,   "NuParentDecX4[4]/D"      );
    fOutputTree->Branch("NuCospibm",      &fNuCospibm,       "NuCospibm/F"       );
    fOutputTree->Branch("NuNorm",         &fNuNorm    ,      "NuNorm/F"                );
    fOutputTree->Branch("NuParentProP4",   fNuParentProP4,   "NuParentProP4[4]/D"      );
    fOutputTree->Branch("NuParentProX4",   fNuParentProX4,   "NuParentProX4[4]/D"      );
    fOutputTree->Branch("NuCospi0bm",     &fNuCospi0bm,      "NuCospi0bm/F"       );
    fOutputTree->Branch("NuParentProNVtx",&fNuParentProNVtx, "NuParentProNVtx/I"       );  //obsolete
    fOutputTree->Branch("NuRnu",          &fNuRnu,           "NuRnu/F"       );
    fOutputTree->Branch("NuXnu",           fNuXnu,           "NuXnu[2]/F"         );
    fOutputTree->Branch("NuIdfd",         &fNuIdfd,          "NuIdfd/I"       );
    fOutputTree->Branch("NuGipart",       &fNuGipart,        "NuGipart/I"              );
    fOutputTree->Branch("NuGpos0",         fNuGpos0,         "NuGpos0[3]/F"            );
    fOutputTree->Branch("NuGvec0",         fNuGvec0,         "NuGvec0[3]/F"            );
    fOutputTree->Branch("NuGamom0",       &fNuGamom0,        "NuGamom0/F"              );


    // jnubeam >= 2010d
    if (2010.35 < fluxVersion) {
      fOutputTree->Branch("NuNg",           &fNuNg,            "NuNg/I"                  );
      fOutputTree->Branch("NuGp",            fNuGp,            "NuGp[NuNg][3]/F"         );
      fOutputTree->Branch("NuGcosbm",        fNuGcosbm,        "NuGcosbm[NuNg]/F"        );
      fOutputTree->Branch("NuGv",            fNuGv,            "NuGv[NuNg][3]/F"         );
      fOutputTree->Branch("NuGpid",          fNuGpid,          "NuGpid[NuNg]/I"          );
      fOutputTree->Branch("NuGmec",          fNuGmec,          "NuGmec[NuNg]/I"          );

      fOutputTree->Branch("NuEnusk",        &fNuEnusk    ,     "NuEnusk/F"               );
      fOutputTree->Branch("NuNormsk",       &fNuNormsk    ,    "NuNormsk/F"              );
      fOutputTree->Branch("NuAnorm",        &fNuAnorm    ,     "NuAnorm/F"               );

      // jnubeam >= 2011a
      if (2011 <= fluxVersion) {
	fOutputTree->Branch("NuGmat",          fNuGmat,          "NuGmat[NuNg]/I"         );
	fOutputTree->Branch("NuGdistc",        fNuGdistc,        "NuGdistc[NuNg]/F"       );
	fOutputTree->Branch("NuGdistal",       fNuGdistal,       "NuGdistal[NuNg]/F"      );
	fOutputTree->Branch("NuGdistti",       fNuGdistti,       "NuGdistti[NuNg]/F"      );
	fOutputTree->Branch("NuGdistfe",       fNuGdistfe,       "NuGdistfe[NuNg]/F"      );
      }

      // Beam parameters
      fOutputTree->Branch("NuVersion",      &fNuVersion,       "NuVersion/F"            );
      fOutputTree->Branch("NuTuneid",       &fNuTuneid,        "NuTuneid/I"             );
      fOutputTree->Branch("NuNtrig",        &fNuNtrig,         "NuNtrig/I"              );
      fOutputTree->Branch("NuPint",         &fNuPint,          "NuPint/I"               );
      fOutputTree->Branch("NuBpos",          fNuBpos,          "NuBpos[2]/F"            );
      fOutputTree->Branch("NuBtilt",         fNuBtilt,         "NuBtilt[2]/F"           );
      fOutputTree->Branch("NuBrms",          fNuBrms,          "NuBrms[2]/F"            );
      fOutputTree->Branch("NuEmit",          fNuEmit,          "NuEmit[2]/F"            );
      fOutputTree->Branch("NuAlpha",         fNuAlpha,         "NuAlpha[2]/F"           );
      fOutputTree->Branch("NuHcur",          fNuHcur,          "NuHcur[3]/F"            );
      fOutputTree->Branch("NuRand",        &fNuRand,           "NuRand/I"               );
      //      fOutputTree->Branch("NuRseed",        fNuRseed,          "NuRseed[2]/I"           );
    }

  } else {
    std::cerr << "Unknown output format " << fFormat << std::endl;
  }

  return;
}

void TNeutOutput::SetTreeWeight(double potToGenerate) {
  
  if (fOutputTree) fOutputTree->SetWeight(potToGenerate);
  else {
    std::cerr << "TNeutOutput::SetTreeWeight() Error: Attempting to set tree weight to un-initialized pointer" << std::endl;
    exit (-1);
  }

}

void TNeutOutput::FillTree(Int_t                eventNumber,
			   TNuVector*           theNuVector,
			   TNuVertex*           theNuVertex,
			   //			   TNeutEventGenerator* theNuFinalState){
			   TNuEventGenerator* theNuFinalState){
  if(fFormat == 0){
    ClearNeutTree(); 
    FillNeutTree(eventNumber, theNuVector, theNuVertex, theNuFinalState);
  } else if(fFormat == 1){
    ClearRooTrackerTree();
    FillRooTrackerTree(eventNumber, theNuVector, theNuVertex, theNuFinalState);
  } else {
    std::cerr << "Unknown output format " << fFormat << std::endl;
  }

}

void TNeutOutput::ClearNeutTree(){
  fNev    = -999;
  fPos[0] = -999;
  fPos[1] = -999;
  fPos[2] = -999;
  fMode   = -999;
  fNumnu  = 0;
  fNpar = 0;
  for(int ip = 0; ip < 100; ip++){
    fAbspnu[ip]  = -999;
    fPnu[ip][0]  = -999;
    fPnu[ip][1]  = -999;
    fPnu[ip][2]  = -999;
    fIpv[ip]     = -999;
    fIcrnv[ip]   = -999;
    fPmomv[ip][0] = -999;
    fPmomv[ip][1] = -999;
    fPmomv[ip][2] = -999;
  }
  
  for(int ip = 0; ip <3; ip++){
    fxpi[ip]   = -999;
    fnpi[ip]   = -999;
    fxpi0[ip]  = -999;
    fnpi0[ip]  = -999;
  }


  fcospibm  = -999;
  fppi      = -999;
  fppid     = -999;
  fcospi0bm = -999;
  fppi0     = -999;
   
  return;
}

void TNeutOutput::ClearRooTrackerTree(){
  fEvtCode  = 0;
  fEvtNum   = -999999;
  fNEibound = -999999;
  fEvtXSec  = -999999;
  fEvtDXSec = -999999;
  fEvtWght  = -999999;
  fEvtProb  = -999999;
  fStdHepN  = 0;
  for(int ip = 0; ip < 4; ip++){
    fEvtVtx[ip] = -999999;
  }

  for(int ip = 0; ip < kNPmax; ip++){
    fStdHepPdg[ip]    = -999999;
    fStdHepStatus[ip] = -999999;

    for(int ip2 = 0; ip2 < 4; ip2++){
      fStdHepX4[ip][ip2] = -999999;
      fStdHepP4[ip][ip2] = -999999;
    }
    for(int ip2 = 0; ip2 < 3; ip2++){
      fStdHepPolz[ip][ip2] = -999999;
    }
    fStdHepFd[ip] = -999999;
    fStdHepLd[ip] = -999999;
    fStdHepFm[ip] = -999999;
    fStdHepLm[ip] = -999999;
  }

  fNEnvc = 0;
  for(int ip = 0; ip < kNEmaxvc; ip++){
    fNEipvc[ip] = -999999;
    fNEiorgvc[ip] = -999999;
    fNEiflgvc[ip] = -999999;
    fNEicrnvc[ip] = -999999;
    for(int ip2 = 0; ip2 < 3; ip2++){
      fNEpvc[ip][ip2] = -999999;
    }
  }

  fNEneutmode = -999999;
  fNEnvert = 0;
  for(int iv = 0; iv < kNEmaxvert; iv++){
    fNEiflgvert[iv] = -999999;
    for(int ip2 = 0; ip2 < 3; ip2++){
      fNEposvert[iv][ip2] = -999999;
    }
  }

  fNEnvcvert = 0;  
  for(int ip = 0; ip < kNEmaxvertp; ip++){
    fNEabspvert[ip] = -999999;
    fNEabstpvert[ip] = -999999;
    fNEipvert[ip] = -999999;
    fNEiverti[ip] = -999999;
    fNEivertf[ip] = -999999;
    for(int ip2 = 0; ip2 < 3; ip2++){
      fNEdirvert[ip][ip2] = -999999;
    }
  }

  fNuParentPdg     = -999999;
  fNuParentDecMode = -999999;
  for(int ip = 0; ip < 4; ip++){
    fNuParentDecP4[ip] = -999999;
    fNuParentDecX4[ip] = -999999;
    fNuParentProP4[ip] = -999999;
    fNuParentProX4[ip] = -999999;
  }

  fNuCospibm = -999999;

  fNuNorm = -999999;

  fNuCospi0bm = -999999;

  fNuParentProNVtx = -999999;

  fNuRnu = -999999;
  for(int ip = 0; ip < 2; ip++){
    fNuXnu[ip] = -999999;
  }      

  fNuIdfd = -999999;

  fNuGipart = -999999;
  for(int ip = 0; ip < 3; ip++){
    fNuGpos0[ip] = -999999;
    fNuGvec0[ip] = -999999;
  }
  fNuGamom0 = -999999;

  fNuNg = 0;
  for(int ip = 0; ip < kNgmax; ip++){
    fNuGpid[ip] = -999999;
    fNuGmec[ip] = -999999;
    for(int ip2 = 0; ip2 < 3; ip2++){
      fNuGv[ip][ip2] = -999999;
      fNuGp[ip][ip2] = -999999;
    }
    fNuGcosbm[ip] = -999999;

    fNuGmat[ip] = -999999;
    fNuGdistc[ip] = -999999;
    fNuGdistal[ip] = -999999;
    fNuGdistti[ip] = -999999;
    fNuGdistfe[ip] = -999999;

  }
  
  fNuEnusk = -999999;
  fNuNormsk = -999999;
  fNuAnorm = -999999;




  fNuVersion = -999999;
  fNuTuneid = -999999;
  fNuNtrig = -999999;
  fNuPint = -999999;
  for(int ip = 0; ip < 2; ip++){
    fNuBpos[ip] = -999999;
    fNuBtilt[ip] = -999999;
    fNuBrms[ip] = -999999;
    fNuEmit[ip] = -999999;
    fNuAlpha[ip] = -999999;
  }
  for(int ip2 = 0; ip2 < 3; ip2++)
    fNuHcur[ip2] = -999999;

  fNuRand = -999999;
  fNuRseed[0] = -999999;
  fNuRseed[1] = -999999;
  
  if(fNuFileName) { 
    delete fNuFileName;
    fNuFileName = NULL;
  }
  fNuFluxEntry = -999999; 
}


void TNeutOutput::FillNeutTree(Int_t                eventNumber,
			       TNuVector*           theNuVector,
			       TNuVertex*           theNuVertex,
			       TNeutEventGenerator* theNuFinalState){

  fNev   = eventNumber;

  fPos[0]   = theNuVertex->getVertex().X()/CLHEP::cm;
  fPos[1]   = theNuVertex->getVertex().Y()/CLHEP::cm;
  fPos[2]   = theNuVertex->getVertex().Z()/CLHEP::cm;

  fMode     = theNuFinalState->getInteractionMode();

  // initial paricle stack
  fNumnu = theNuFinalState->getNNEParticles();
  for(int ip = 0; ip < fNumnu; ip++){
    fIpnu[ip]   = theNuFinalState->getNEID(ip);
    Float_t* p = theNuFinalState->getNEMomentum(ip);
    fAbspnu[ip] = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])/CLHEP::GeV;
    fPnu[ip][0] = p[0]/CLHEP::GeV;
    fPnu[ip][1] = p[1]/CLHEP::GeV;
    fPnu[ip][2] = p[2]/CLHEP::GeV;
  }

  // particles to track
  fNpar = theNuFinalState->getNVCParticles();
  Int_t ic = 0;
  for(int ip = 0; ip < fNpar; ip++){

    // this bit of code leaves ip=1 open for the target nucleus
    ic = ip;
    if(ip > 0)ic = ip + 1;

    fIpv[ic]   = theNuFinalState->getVCID(ip); 
    fIcrnv[ic] = theNuFinalState->getVCChase(ip);

    Float_t* p = theNuFinalState->getVCMomentum(ip);
    fPmomv[ic][0] = p[0]/CLHEP::MeV;   
    fPmomv[ic][1] = p[1]/CLHEP::MeV;   
    fPmomv[ic][2] = p[2]/CLHEP::MeV;   
  }

  // throw in the target nucleus here
  fIpv[1]   = theNuVertex->getPDGofNucleus();  
  fIcrnv[1] = 0;
  
  // nucleus is at rest;
  fPmomv[1][0] = 0; 
  fPmomv[1][1] = 0;
  fPmomv[1][2] = 0;

  fNpar++;

  // now to the neutrino vector information from jnubeam
  TVector3 theDecayPos      = theNuVector->getDecayPos();
  TVector3 theDecayMomentum = theNuVector->getDecayMomentum();

  fxpi[0] =  theDecayPos.X()/CLHEP::cm;     
  fxpi[1] =  theDecayPos.Y()/CLHEP::cm;    
  fxpi[2] =  theDecayPos.Z()/CLHEP::cm; 

  fnpi[0]  = theDecayMomentum.Unit().X();
  fnpi[1]  = theDecayMomentum.Unit().Y();
  fnpi[2]  = theDecayMomentum.Unit().Z();

  fcospibm = theNuVector->getCosToBeamAtDecay();
  fppi     = theDecayMomentum.Mag()/CLHEP::GeV;    
  fppid    = theNuVector->getDecayParticleID();   

  TVector3 theProductionMomentum = theNuVector->getProductionMomentum();
  TVector3 theProductionPosition = theNuVector->getProductionPos();

  fxpi0[0] = theProductionPosition.X()/CLHEP::cm;     
  fxpi0[1] = theProductionPosition.Y()/CLHEP::cm;     
  fxpi0[2] = theProductionPosition.Z()/CLHEP::cm;     

  fnpi0[0] = theProductionMomentum.Unit().X();
  fnpi0[1] = theProductionMomentum.Unit().Y();
  fnpi0[2] = theProductionMomentum.Unit().Z();

  fcospi0bm = theNuVector->getCosToBeamAtProduction();
  fppi0     = theProductionMomentum.Mag()/CLHEP::GeV;

  fOutputTree->Fill();

}

void TNeutOutput::FillRooTrackerTree(Int_t                eventNumber,
				     TNuVector*           theNuVector,
				     TNuVertex*           theNuVertex,
				     //TNeutEventGenerator* theNuFinalState){
				     TNuEventGenerator* theNuFinalState){

  /// The generator-specific string with the 'event code'
  fEvtCode = new TObjString(Form("%d", theNuFinalState->getInteractionMode()) );
  
  /// The sequence number of the event (the event number).
  fEvtNum = eventNumber;
  
  /// Bound target nucleon or not
  fNEibound = theNuFinalState->getIbound();
  
  /// The cross section for the event (1E-38 cm2)
  fEvtXSec = theNuFinalState->getTotcrs();
  
  /// The differential cross section for the event kinematics 
  /// (1E-38 cm2/{K^n})
  //fEvtDXSec = -1; // not filled currently
  
  /// The weight for the event
  //fEvtWght  = 1;    // events are not weighted
  
  /// The probability for the event (given the cross section, path lengths,
  /// etc.).
  //fEvtProb = -1;   // not filled currently
  
  /// The event vertex position in detector coord syst (in meters and
  /// seconds).
  fEvtVtx[0]   = theNuVertex->getVertex().X()/CLHEP::m;
  fEvtVtx[1]   = theNuVertex->getVertex().Y()/CLHEP::m;
  fEvtVtx[2]   = theNuVertex->getVertex().Z()/CLHEP::m;
  fEvtVtx[3]   = theNuVertex->getVertex().T()/CLHEP::s;
  
  Int_t ic = 0;
  for(int ip = 0; ip < theNuFinalState->getNVCParticles(); ip++){

    // ugly counter gymnastics to sneak in the nucleus information
    // in position 1 following neutrino information. Thus we need
    // to leave a gap at 1
    if (ip==1) ic++;

    fStdHepPdg[ic]    = theNuFinalState->getVCID(ip);
    
    Int_t neutFlag = theNuFinalState->getVCFlags(ip);
    
    // Attempt to convert neut flags according to these conventions
    // The GENIE generator approximately follows the HEPEVT status codes.
    // As of July 2008, the status values found the GENIE source code are:
    //   - -1 -- Undefined particle
    //   -  0 -- An initial state particle.
    //   -  1 -- A stable final state particle to be tracked.
    //   -  2 -- An intermediate particle that should not be tracked.
    //   -  3 -- A particle which decayed and should not be tracked.  If 
    //            this particle produced daughters to be tracked, they will
    //            have a state of 1.

    if(     neutFlag   == -1){// neut initial state
      fStdHepStatus[ic] =  0; // HEPEVT initial state particle
    } 
    else if(neutFlag   ==  0){// neut  "determined later" =  final state
      fStdHepStatus[ic] =  1; // HEPEVT stable final state particle
    }
    else if(neutFlag   ==  2){// neut "escape from detector" (e.g. neutrino)
      fStdHepStatus[ic] =  1; // HEPEVT stable final state particle
    } else {
      continue;
    }

    // Rest is deprecated (Full info in NEUT vcwork common block below)
    /*
    else if(neutFlag   ==  1){// neut "decay to other particle"
      fStdHepStatus[ic] =  3; // HEPEVT intermediate particle which decayed
    }
    else if(neutFlag   ==  3){// neut absorption (may produce secondaries)
      fStdHepStatus[ic] =  3; // HEPEVT intermediate particle which decayed
    }
    else if(neutFlag   ==  4){// neut charge exchange (may produce secondaries)
      fStdHepStatus[ic] =  3; // HEPEVT intermediate particle which decayed
    }
    else if(neutFlag   ==  5){// neut "stop and ignore in MC"
      fStdHepStatus[ic] =  2; // HEPEVT intermediate that should not be tracked
    }
    else if(neutFlag   ==  6){// neut "EM shower" ???
      fStdHepStatus[ic] =  3; // HEPEVT intermediate particle which decayed
    }
    else if(neutFlag   ==  7){// neut particle production
      fStdHepStatus[ic] =  3; // HEPEVT intermediate particle which decayed
    }
    else if(neutFlag   ==  8){// neut inelastic scatter
      fStdHepStatus[ic] =  3; // HEPEVT intermediate particle which decayed
    } else {
      std::cerr << "Unknown neut status code " << neutFlag << std::endl;
    }
    */

    /// The position (x, y, z, t) (fm, second) of the particle in the nuclear
    /// frame
    //fStdHepX4[ic][0] = -999999; // not filled in Neut currently
    //fStdHepX4[ic][1] = -999999; // not filled in Neut currently
    //fStdHepX4[ic][2] = -999999; // not filled in Neut currently
    //fStdHepX4[ic][3] = -999999; // not filled in Neut currently

    /// The 4-momentum (px, py, pz, E) of the particle in the LAB frame (GeV)
    Float_t* p  = theNuFinalState->getVCMomentum(ip);
    Float_t   m = 0;
    if(TParticlePDG* theParticle = fPDGDatabase.GetParticle(fStdHepPdg[ic]) ){
      m = theParticle->Mass()*CLHEP::GeV; //TParticlePDG returns masses in GeV
    } else {
      //std::cerr << "Unknown particle: PDG " << fStdHepPdg[ic] 
      //<< " Setting mass to zero " << std::endl;
      m = 0;
    }
    
    // momentum returned in MeV
    fStdHepP4[ic][0] = p[0]/CLHEP::GeV;
    fStdHepP4[ic][1] = p[1]/CLHEP::GeV;
    fStdHepP4[ic][2] = p[2]/CLHEP::GeV;
    fStdHepP4[ic][3] = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + m*m)/CLHEP::GeV;
      
    /// The particle polarization vector.
    // polarizations not reported by Neut
    //fStdHepPolz  [ic][0] = -999999;
    //fStdHepPolz  [ic][1] = -999999;
    //fStdHepPolz  [ic][2] = -999999;
      
    // determine first and last daughter particles of this particle 
    // DEPRECATED (will be included in full NEUT vcwork common block later)
    /*
    Int_t firstDaughter = -1;
    Int_t lastDaughter  = -1;
       
    for(int ip2 = 0; ip2 < theNuFinalState->getNVCParticles(); ip2++){
      Int_t mother = theNuFinalState->getVCOrigin(ip2);
      if(mother == ip && firstDaughter == -1){
	firstDaughter = ip2;
      }
      if(mother == ip && ip2 > lastDaughter){
	lastDaughter = ip2;
      }
    }
 
    /// The index of the first daughter of the particle in the arrays.
    fStdHepFd[ic] = firstDaughter;
    fStdHepLd[ic] = lastDaughter;
    
    // neut particles only have one mother
    fStdHepFm[ic] = theNuFinalState->getVCOrigin(ip);
    fStdHepLm[ic] = theNuFinalState->getVCOrigin(ip);
    */
   
    //fStdHepFd[ic] = -999999;
    //fStdHepLd[ic] = -999999;
    //fStdHepFm[ic] = -999999;
    //fStdHepLm[ic] = -999999;

    ic++;
    fStdHepN++;
    

  } // end of loop over particle stack
  
  
  //
  // Sneak in the nucleus information on the stack
  //
  fStdHepPdg[1]    = theNuVertex->getPDGofNucleus();  
  fStdHepStatus[1] = 0; // HEPEVT initial state particle  
  
  fStdHepX4[1][0]  = 0; // not filled in Neut currently
  fStdHepX4[1][1]  = 0; // not filled in Neut currently
  fStdHepX4[1][2]  = 0; // not filled in Neut currently
  fStdHepX4[1][3]  = 0; // not filled in Neut currently

  fStdHepP4[1][0] = 0;
  fStdHepP4[1][1] = 0;
  fStdHepP4[1][2] = 0;
  fStdHepP4[1][3] = theNuVertex->getA();
  
  /// The particle polarization vector.
  // polarizations not reported by Neut
  fStdHepPolz[1][0] = 0;
  fStdHepPolz[1][1] = 0;
  fStdHepPolz[1][2] = 0;  

  // do not set mother/daughter information for nucleus
  // note that this is different from Genie (which tracks nucleons
  // as daughters of the nucleus
  fStdHepFd[1] = -1;
  fStdHepLd[1] = -1;
  fStdHepFm[1] = -1;
  fStdHepLm[1] = -1;
  fStdHepN++;


  // Fill some NEUT kinematic information
  fNEcrsx  = theNuFinalState->getCrsX();
  fNEcrsy  = theNuFinalState->getCrsY();
  fNEcrsz  = theNuFinalState->getCrsZ();
  fNEcrsphi= theNuFinalState->getCrsPhi(); 

  // Fill native NEUT VCWORK information
  fNEneutmode = theNuFinalState->getInteractionMode();
  fNEnvc = theNuFinalState->getNVCParticles();
  for(int ip = 0; ip < fNEnvc; ip++){
    fNEipvc[ip] = theNuFinalState->getVCID(ip);
    fNEiorgvc[ip] = theNuFinalState->getVCOrigin(ip);
    fNEiflgvc[ip] = theNuFinalState->getVCFlags(ip);
    fNEicrnvc[ip] = theNuFinalState->getVCChase(ip);

    Float_t* p  = theNuFinalState->getVCMomentum(ip);
    for(int ip2 = 0; ip2 < 3; ip2++){
      fNEpvc[ip][ip2] = p[ip2];  
    }
  }


  // Fill NEUT FSIHIST information
  fNEnvert = theNuFinalState->getNFSVertices();
  for(int iv = 0; iv < fNEnvert; iv++){
    fNEiflgvert[iv] = theNuFinalState->getFSIntFlag(iv);
  
    Float_t* pos  = theNuFinalState->getFSPos(iv);
    for(int ip2 = 0; ip2 < 3; ip2++){
      fNEposvert[iv][ip2] = pos[ip2];  
    }
  }  

  
  fNEnvcvert = theNuFinalState->getNFSParticles();
  for(int ip = 0; ip < fNEnvcvert; ip++){
    fNEabspvert[ip] = theNuFinalState->getFSMomentumLab(ip);
    fNEabstpvert[ip] = theNuFinalState->getFSMomentumNuc(ip);
    fNEipvert[ip] = theNuFinalState->getFSID(ip);
    fNEiverti[ip] = theNuFinalState->getFSVertexInit(ip);
    fNEivertf[ip] = theNuFinalState->getFSVertexFin(ip);

    Float_t* dir  = theNuFinalState->getFSDir(ip);
    for(int ip2 = 0; ip2 < 3; ip2++){
      fNEdirvert[ip][ip2] = dir[ip2];
    }
  }

  // Fill in the Beam Information
  // beam information

  if(fNuFileName) { 
    delete fNuFileName;
    fNuFileName = 0;
  }
  fNuFileName = new TObjString((theNuVector->getNuFileName()).Data());
  fNuFluxEntry = theNuVector->getNuFluxEntry();


  // The PDG code of the particle which created the parent neutrino.
  fNuParentPdg = fPDGDatabase.ConvertGeant3ToPdg(theNuVector->getDecayParticleID());
  
  /// The interaction mode at the vertex which created the parent neutrino.
  /// This is normally the decay mode of the parent particle.
  fNuParentDecMode =  theNuVector->getDecayMode();
  
  /// The 4 momentum of the particle at the vertex which created the parent
  /// neutrino.  This is normally the momentum of the parent particle at the
  /// decay point.
  TVector3 theDecayMomentum = theNuVector->getDecayMomentum();
  // mass from PDG database is returned in GeV/c^2.
  Double_t mParent          = (fPDGDatabase.GetParticle(fNuParentPdg))->Mass();

  fNuParentDecP4[0] = theDecayMomentum.X()/CLHEP::GeV;
  fNuParentDecP4[1] = theDecayMomentum.Y()/CLHEP::GeV;
  fNuParentDecP4[2] = theDecayMomentum.Z()/CLHEP::GeV;
  fNuParentDecP4[3] = sqrt(fNuParentDecP4[0] * fNuParentDecP4[0] +
			   fNuParentDecP4[1] * fNuParentDecP4[1] +
			   fNuParentDecP4[2] * fNuParentDecP4[2] +
			   mParent * mParent);
  
  /// The position (meters, seconds) of the vertex at which the neutrino was
  /// created.  This uses the target as the origin.
  TVector3 theDecayPos = theNuVector->getDecayPos();
  fNuParentDecX4[0] = theDecayPos.X()/CLHEP::cm;
  fNuParentDecX4[1] = theDecayPos.Y()/CLHEP::cm;
  fNuParentDecX4[2] = theDecayPos.Z()/CLHEP::cm;
  fNuParentDecX4[3] = 0/CLHEP::s;
    
  fNuCospibm = theNuVector->getCosToBeamAtDecay();

  /// The momentum (GeV) of the parent particle at it's production point.
  TVector3 theProductionMomentum = theNuVector->getProductionMomentum();
  
  fNuParentProP4[0] = theProductionMomentum.X()/CLHEP::GeV;
  fNuParentProP4[1] = theProductionMomentum.Y()/CLHEP::GeV;
  fNuParentProP4[2] = theProductionMomentum.Z()/CLHEP::GeV;
  fNuParentProP4[3] = sqrt(fNuParentProP4[0] * fNuParentProP4[0] +
			   fNuParentProP4[1] * fNuParentProP4[1] +
			   fNuParentProP4[2] * fNuParentProP4[2] +
			   mParent * mParent);
  
  /// The position (meters, seconds) of the parent particle at it's
  /// production point.  This uses the target as the origin.
  TVector3 theProdPos = theNuVector->getProductionPos();
  fNuParentProX4[0] = theProdPos.X()/CLHEP::cm;
  fNuParentProX4[1] = theProdPos.Y()/CLHEP::cm;
  fNuParentProX4[2] = theProdPos.Z()/CLHEP::cm;
  fNuParentProX4[3] = 0/CLHEP::s;   // currently not filled by jnubeam

  fNuCospi0bm = theNuVector->getCosToBeamAtProduction();

  /// The vertex ID of the parent particle vertex. (obsolete)
  fNuParentProNVtx  = theNuVector->getNvtx0();

  fNuIdfd  = theNuVector->getIdfd();
    
  fNuRnu = theNuVector->getNuRnu()/CLHEP::cm;
  
  TVector3 theNuXnu = theNuVector->getNuPos();
  fNuXnu[0] = theNuXnu.X()/CLHEP::cm;
  fNuXnu[1] = theNuXnu.Y()/CLHEP::cm;
  
  fNuGipart = theNuVector->getPrimaryParticleID();

  // primary particle starting point
  TVector3 thePrimaryParticlePos = theNuVector->getPrimaryParticlePos();
  fNuGpos0[0] =  thePrimaryParticlePos.X()/CLHEP::cm;
  fNuGpos0[1] =  thePrimaryParticlePos.Y()/CLHEP::cm;
  fNuGpos0[2] =  thePrimaryParticlePos.Z()/CLHEP::cm;


  // primary particle direction at starting point
  TVector3 thePrimaryParticleMomentum = theNuVector->getPrimaryParticleMomentum();
  fNuGvec0[0] = thePrimaryParticleMomentum.Unit().X();
  fNuGvec0[1] = thePrimaryParticleMomentum.Unit().Y();
  fNuGvec0[2] = thePrimaryParticleMomentum.Unit().Z();

  // momentum of the primary particle at the starting point
  fNuGamom0 = thePrimaryParticleMomentum.Mag()/CLHEP::GeV;

  fNuNorm = theNuVector->getNuWeight();

  fNuVersion = theNuVector->getNuVersion(); 

  // jnubeam >= 2010d
  if (2010.35<fNuVersion) {

    // Number of interaction steps 
    fNuNg = theNuVector->getNumIntSteps();

    for(int ip = 0; ip < fNuNg; ip++){
     
      TVector3 theIntermediateParticleMomentum = theNuVector->getIntermediateParticleMomentum(ip);
      fNuGp[ip][0] = theIntermediateParticleMomentum.X()/CLHEP::GeV;
      fNuGp[ip][1] = theIntermediateParticleMomentum.Y()/CLHEP::GeV;
      fNuGp[ip][2] = theIntermediateParticleMomentum.Z()/CLHEP::GeV;

      fNuGcosbm[ip] = theNuVector->getCosToBeamOfIntermediate(ip); 
  
      TVector3 theIntermediateParticlePos = theNuVector->getIntermediateParticlePos(ip);
      fNuGv[ip][0] = theIntermediateParticlePos.X()/CLHEP::cm;
      fNuGv[ip][1] = theIntermediateParticlePos.Y()/CLHEP::cm;
      fNuGv[ip][2] = theIntermediateParticlePos.Z()/CLHEP::cm;

      // Intermediate particle ID
      fNuGpid[ip] = theNuVector->getIntermediateParticleID(ip);
  
      // Intermediate particle Interaction mechanism
      fNuGmec[ip] = theNuVector->getIntermediateParticleMec(ip); 

      // Out of target information
      // jnubeam >= 2011a
      if (2011 <= fNuVersion) {
	fNuGmat[ip] = theNuVector->getNuGmat(ip); 
	fNuGdistc[ip] = theNuVector->getNuGdistc(ip)/CLHEP::cm; 
	fNuGdistal[ip] = theNuVector->getNuGdistal(ip)/CLHEP::cm; 
	fNuGdistti[ip] = theNuVector->getNuGdistti(ip)/CLHEP::cm; 
	fNuGdistfe[ip] = theNuVector->getNuGdistfe(ip)/CLHEP::cm; 
      }
    }
 
    fNuEnusk = theNuVector->getNuEnergySK()/CLHEP::GeV;

    fNuNormsk = theNuVector->getNuWeightSK();
    fNuAnorm = theNuVector->getNuWeightAccept(); 
 


    fNuTuneid = theNuVector->getNuTuneid();

    fNuNtrig = theNuVector->getNuNtrig();
  
    fNuPint = theNuVector->getNuPint();

    for (Int_t i=0; i<2; i++) {
      fNuBpos[i] = theNuVector->getNuBpos(i);
      fNuBtilt[i] = theNuVector->getNuBtilt(i);
      
      fNuBrms[i] = theNuVector->getNuBrms(i);
      
      fNuEmit[i] = theNuVector->getNuEmit(i);
      fNuAlpha[i] = theNuVector->getNuAlpha(i);
    }
    for (Int_t i=0; i<3; i++) {
      fNuHcur[i] = theNuVector->getNuHcur(i);
    }

    fNuRand = theNuVector->getNuRand();
    for (Int_t i=0; i<2; i++) 
      fNuRseed[i] = theNuVector->getNuRseed(i);
  }


  fOutputTree->Fill();
  
}
