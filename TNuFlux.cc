#include "TNuFlux.h"

TNuFlux::TNuFlux(){
  initialize();
}

TNuFlux::~TNuFlux(){
}

TNuFlux::TNuFlux(TString fileName, TString treeKey, TString verKey, TString ndKey){
  setTreeKey(treeKey);
  setVerTreeKey(verKey); 
  setNDTreeKey(ndKey); 

  fFileName = fileName;

  // Pre-initialization to check version of first file
  fFileIndex = 0;
  nHaddFluxFiles = 0;
  version = tuneid = ftuneid = 999;
  checkVerTree(fileName);
  
  initialize();

}


void TNuFlux::initialize() {

  nChainFluxFiles = 0;
  fDetector = 0;
  fNuType   = 0;
  fVectorIndex = 0;
  fVectorTotal = 0;
  fRewind      = false;

  // >2010a
  fNuFluxChain.SetBranchAddress("Enu",     &fEnu,      &b_Enu     );
  fNuFluxChain.SetBranchAddress("ppid",    &fppid,     &b_ppid    );
  fNuFluxChain.SetBranchAddress("mode",    &fmode,     &b_mode    );
  fNuFluxChain.SetBranchAddress("ppi",     &fppi,      &b_ppi     );
  fNuFluxChain.SetBranchAddress("xpi",      fxpi,      &b_xpi     );
  fNuFluxChain.SetBranchAddress("npi",      fnpi,      &b_npi     );
  fNuFluxChain.SetBranchAddress("cospibm", &fcospibm,  &b_cospibm );
  fNuFluxChain.SetBranchAddress("norm",    &fnorm,     &b_norm    );
  fNuFluxChain.SetBranchAddress("nvtx0",   &fnvtx0,    &b_nvtx0   );
  fNuFluxChain.SetBranchAddress("ppi0",    &fppi0,     &b_ppi0    );
  fNuFluxChain.SetBranchAddress("xpi0",     fxpi0,     &b_xpi0    );
  fNuFluxChain.SetBranchAddress("npi0",     fnpi0,     &b_npi0    );
  fNuFluxChain.SetBranchAddress("cospi0bm",&fcospi0bm, &b_cospi0bm);
  fNuFluxChain.SetBranchAddress("rnu",     &frnu,      &b_rnu     );
  fNuFluxChain.SetBranchAddress("xnu",     &fxnu,      &b_xnu     );
  fNuFluxChain.SetBranchAddress("ynu",     &fynu,      &b_ynu     );
  fNuFluxChain.SetBranchAddress("nnu",      fnnu,      &b_nnu     );
  fNuFluxChain.SetBranchAddress("idfd",    &fidfd,     &b_idfd    );

  fNuFluxChain.SetBranchAddress("gipart",  &fgipart,   &b_gipart  );
  fNuFluxChain.SetBranchAddress("gpos0" ,  fgpos0,     &b_gpos0   );
  fNuFluxChain.SetBranchAddress("gvec0" ,  fgvec0,     &b_gvec0   );
  fNuFluxChain.SetBranchAddress("gamom0",  &fgamom0,   &b_gamom0  );

  // jnubeam > 2010d
  if (2010.35 < version) {
    fNuFluxChain.SetBranchAddress("ng", &fng, &b_ng);
    fNuFluxChain.SetBranchAddress("gpx", fgpx, &b_gpx);
    fNuFluxChain.SetBranchAddress("gpy", fgpy, &b_gpy);
    fNuFluxChain.SetBranchAddress("gpz", fgpz, &b_gpz);
    fNuFluxChain.SetBranchAddress("gcosbm", fgcosbm, &b_gcosbm);
    fNuFluxChain.SetBranchAddress("gvx", fgvx, &b_gvx);
    fNuFluxChain.SetBranchAddress("gvy", fgvy, &b_gvy);
    fNuFluxChain.SetBranchAddress("gvz", fgvz, &b_gvz);
    fNuFluxChain.SetBranchAddress("gpid", fgpid, &b_gpid);
    fNuFluxChain.SetBranchAddress("gmec", fgmec, &b_gmec);
    fNuFluxChain.SetBranchAddress("Enusk", &fEnusk, &b_Enusk);
    fNuFluxChain.SetBranchAddress("normsk", &fnormsk, &b_normsk);
    fNuFluxChain.SetBranchAddress("anorm", &fanorm, &b_anorm);

    //jnubeam >= 2011a
    if (2011 <= version) {
      fNuFluxChain.SetBranchAddress("gmat", fgmat, &b_gmat);
      fNuFluxChain.SetBranchAddress("gdistc", fgdistc, &b_gdistc);
      fNuFluxChain.SetBranchAddress("gdistal", fgdistal, &b_gdistal);
      fNuFluxChain.SetBranchAddress("gdistti", fgdistti, &b_gdistti);
      fNuFluxChain.SetBranchAddress("gdistfe", fgdistfe, &b_gdistfe);
    }

    // Version tree
    fVerTree.SetBranchAddress("version", &fversion, &b_version);
    fVerTree.SetBranchAddress("tuneid", &ftuneid, &b_tuneid);
    fVerTree.SetBranchAddress("ntrig", &fntrig, &b_ntrig);
    fVerTree.SetBranchAddress("pint", &fpint, &b_pint);
    fVerTree.SetBranchAddress("bpos", fbpos, &b_bpos);
    fVerTree.SetBranchAddress("btilt", fbtilt, &b_btilt);
    fVerTree.SetBranchAddress("brms", fbrms, &b_brms);
    fVerTree.SetBranchAddress("emit", femit, &b_emit);
    fVerTree.SetBranchAddress("alpha", falpha, &b_alpha);
    fVerTree.SetBranchAddress("hcur", fhcur, &b_hcur);
    fVerTree.SetBranchAddress("rand", &frand, &b_rand);
    fVerTree.SetBranchAddress("rseed", frseed, &b_rseed);
  }

  return;
}


void TNuFlux::addToChain(TString theFile){
  fNuFluxChain.Add(Form("%s/%s", theFile.Data(), fKey.Data()));
  fVectorTotal = (Int_t)fNuFluxChain.GetEntries();
  nChainFluxFiles++;
}

int TNuFlux::checkVerTree(TString theFile){

  // Warning: if first file is version <2010d, then this routine does not check 
  // chaining of additional files (this is to suppress error message)
  if (version >= 0) {
    fVerTree.Add(Form("%s/%s", theFile.Data(), fVerKey.Data()));
  }

  fFileIndex++;

  // Store number of files in first file
  static Int_t nHaddFluxFilesFirst = 0;
  
  Int_t nEntries = fVerTree.GetEntries();
  nHaddFluxFiles = nEntries - nHaddFluxFilesFirst;

  // If the tree exists (ver > 2010d), # of entries should increment
  if (nEntries < fFileIndex) {
    // Set negative to skip whole routine on next call
    fversion = -1;
    ftuneid = -1;
    
    // Things done at initialization only
    if (fFileIndex==1) {
      std::cout << std::endl << "Loading jnubeam flux version < 2010d" << std::endl;
      std::cout << "WARNING: No version control. Make sure all chained flux files are of the same version" << std::endl;
      version = fversion;
      tuneid = ftuneid;
    }
  }
  // Tree exists (ver >= 2010d) so get the variables
  else {

    
    // Things done at initialization only
    if (fFileIndex==1) {
      fVerTree.SetBranchAddress("version", &fversion, &b_version);
      fVerTree.SetBranchAddress("tuneid", &ftuneid, &b_tuneid);

      // Do not start counting yet
      nHaddFluxFilesFirst = nEntries;
    }

    // Check that header of all hadd'ed files is the same
    fVerTree.GetEntry(0);
    Float_t version_prev = fversion;
    Int_t tuneid_prev = ftuneid;
    for (int i=1; i<nEntries; i++) {
      fVerTree.GetEntry(i);
      if (fversion!=version_prev) {
	std::cout << "Version mismatch in hadd'ed flux file #" << i << ": " << fversion << " != " << version_prev << std::endl;
	return -1;
      }
      if (ftuneid!=tuneid_prev) {
	std::cout << "Tuneid mismatch in hadd'ed flux file #" << i << ": " << ftuneid << " != " << tuneid_prev << std::endl;
	return -1;
      }
    }
    
    // Set version of this flux
    if (fFileIndex==1) {
      fVerTree.GetEntry(0);
      version = fversion;
      tuneid = ftuneid;
      
      std::cout << std::endl << "Loading jnubeam flux version = " << version << ", tuneid = " << tuneid << std::endl;

      // Reset the chain so that it has the right number of entries in hadd'ed files
      //fVerTree.Reset();
    }
  }

  // Check version of files being chained together
  if (fFileIndex>1) {
    if (version != fversion) {
      std::cerr << "TNuFlux error: Adding flux file version " << fversion << " to chain of version " << version  << std::endl;
      return -1;
    }
    if (tuneid != ftuneid) {
      std::cerr << "TNuFlux error: Adding flux file tuneid " << ftuneid << " to chain of tuneid " << tuneid << std::endl;
      return -1;
    }
  }

  return 0;
}

void TNuFlux::setTreeKey(TString s){
  fKey = s;
};

TString TNuFlux::getTreeKey(){
  return fKey;
}

void TNuFlux::setVerTreeKey(TString s){
  fVerKey = s;
};

TString TNuFlux::getVerTreeKey(){
  return fVerKey;
}

void TNuFlux::setNDTreeKey(TString s){
  fNDKey = s;
};

TString TNuFlux::getNDTreeKey(){
  return fNDKey;
}

// select specific detector
// 0 = all, 1 = 2 km,  2 = on-axis for 2.5 deg, 3 = INGRID x, 
// 4 = INGRID y (set to 3 for now), 5 = nd280 basket, 6 = nd280 magnet
void TNuFlux::setDetector(Int_t i){
  fDetector = i;

  setPlanePos();
}

Int_t TNuFlux::getDetector(){
  return fDetector;
}

// Grab JNuBeam flux plane coordinates from tree
void TNuFlux::setPlanePos(){
  TChain ndTree;
  ndTree.Add(Form("%s/%s", fFileName.Data(), fNDKey.Data()));

  // Declaration of leaf types
  Int_t           Nfd;
  Float_t         Bxfd[13];   //[Nfd]
  Float_t         Byfd[13];   //[Nfd]
  Float_t         Xfd[13];   //[Nfd]
  Float_t         Yfd[13];   //[Nfd]
  Float_t         Zfd[13];   //[Nfd]
  Float_t         Hfd[13];   //[Nfd]
  Float_t         Vfd[13];   //[Nfd]

  // List of branches
  TBranch        *b_Nfd;   //!
  TBranch        *b_Bxfd;   //!
  TBranch        *b_Byfd;   //!
  TBranch        *b_Xfd;   //!
  TBranch        *b_Yfd;   //!
  TBranch        *b_Zfd;   //!
  TBranch        *b_Hfd;   //!
  TBranch        *b_Vfd;   //!

  ndTree.SetBranchAddress("Nfd", &Nfd, &b_Nfd);
  ndTree.SetBranchAddress("Bxfd", Bxfd, &b_Bxfd);
  ndTree.SetBranchAddress("Byfd", Byfd, &b_Byfd);
  ndTree.SetBranchAddress("Xfd", Xfd, &b_Xfd);
  ndTree.SetBranchAddress("Yfd", Yfd, &b_Yfd);
  ndTree.SetBranchAddress("Zfd", Zfd, &b_Zfd);
  ndTree.SetBranchAddress("Hfd", Hfd, &b_Hfd);
  ndTree.SetBranchAddress("Vfd", Vfd, &b_Vfd);

  ndTree.GetEntry(0);

  // Get z-offset from ND280 origin (warning hardcoded values, corresponding to jnubeam/nubin/nubeam.card)
  for (int ifd=0; ifd<Nfd; ifd++) {
    if (Zfd[ifd]>0) { // Check that z-position makes sense
      fFluxPlanePos[ifd][0] = Xfd[ifd] - (-322.2);
      fFluxPlanePos[ifd][1] = Yfd[ifd] - (-814.6);
      fFluxPlanePos[ifd][2] = Zfd[ifd] - (28010);
    } else { // Negative z-position means flux file does not contain given ND plane
             // Set far downstream so that no interactions are generated in case this
             // plane is accidentally used.
      for (int i=0; i<3; i++)
	fFluxPlanePos[ifd][i] = 999999;
    }
  }  

  // Check if plane is centered at ND280 origin (basket center)
  double epsilon=5e-5;
  
  bool *isCentered = new bool[Nfd];

  std::cout << std::endl;
  for (int ifd=0; ifd<Nfd; ifd++) {

    isCentered[ifd] = 1;
  
    for (int i=0; i<3; i++) {
      if (TMath::Abs(fFluxPlanePos[ifd][i]) > epsilon) isCentered[ifd]=0;
    }
  
    if (!isCentered[ifd]) {
      if (Zfd[ifd]>0) {
	std::cout << "Using ND" << ifd+1 << " flux plane, not centered at ND280 magnet center:" << std::endl;
	std::cout << "(x, y, z) = ";
	for (int i=0; i<3; i++) 
	  std::cout << fFluxPlanePos[ifd][i] << " ";
	std::cout << std::endl;
      }
    } else {
      if (Zfd[ifd]>0)
	std::cout << "Using ND" << ifd+1 << " flux plane at ND280 magnet center" << std::endl;
    }
  }
  std::cout << std::endl;
  delete isCentered;

}


// select neutrino type
// 0 = all , 1 = numu, 2 = numub, 3 = nue, 4 = nueb
void TNuFlux::setNuType(Int_t i){
  fNuType = i;
}

Int_t TNuFlux::getNuType(){
  return fNuType;
}

  
TNuVector* TNuFlux::getVector(){
  
  Long64_t  nb = 0;

  // check if we still have vectors
  if(fVectorIndex >= fVectorTotal){
    std::cout << "TNuFlux::getVector(): Exhausted supply of neutrino vectors on vector number:" 
	      << "\t" 
	      << fVectorIndex << std::endl;
    fVectorIndex = 0;

    if(fRewind){
      rewindTree();
      std::cout << "Rewinding jnubeam tree" << std::endl;
    } else {
      return 0;
    }

  }

  TNuVector* theNuVector = new TNuVector();
  
  do{

    if(fVectorIndex >= fVectorTotal){
      std::cout << "TNuFlux::getVector(): Exhausted supply of neutrino vectors on vector number:" 
		<< "\t" 
		<< fVectorIndex << std::endl;

      if(fRewind){
	std::cout << "Rewinding jnubeam tree " << std::endl;
	rewindTree();
      }	else {
	return 0;
      }
    }

    nb = fNuFluxChain.GetEntry(fVectorIndex);
    
    // Entry in the current (single) flux file
    fNuFluxEntry = fVectorIndex-(int)fNuFluxChain.GetTreeOffset()[(int)fNuFluxChain.GetTreeNumber()];
    
    fVectorIndex++;

        
  } while(!ValidVector(fidfd, fmode/10));

  // Get Current File Name and append tree name
  std::string filename = (fNuFluxChain.GetFile())->GetName();
  std::string::size_type start_pos = filename.rfind("/");
  if (start_pos == std::string::npos) start_pos = 0; else ++start_pos;
  std::string basename = filename.substr(start_pos);
  fCurrentFileName = Form("%s:%s",basename.c_str(),fKey.Data());
  
  FillVector(theNuVector);

  return theNuVector;
}

Bool_t TNuFlux::ValidVector(Int_t detectorID, Int_t neutrinoID){

  // detector id must match  (detector = detectorID);
  // and either we are not selecting neutrino flavor (nuType = 0), or the flavor 
  // matches (nuType = neutrinoID)

  // INGRID ntuple contains both detectorID = 3,4,8,9, but genev and 
  // event_rate accept -d 3 only for now.
  if (detectorID==4 || detectorID==8 || detectorID==9) detectorID=3;

  return (fDetector == detectorID) &&
         (fNuType == 0 || (fNuType == neutrinoID) );
}

void TNuFlux::FillVector(TNuVector* theVector){

  theVector->setNuFileName(fCurrentFileName);
  
  theVector->setEntry(fVectorIndex-1);
  theVector->setNuFluxEntry(fNuFluxEntry);


  theVector->setNuDir(fnnu[0], fnnu[1], fnnu[2]); 

  theVector->setNuEnergy(fEnu * CLHEP::GeV);    

  Int_t neutrinoID = -1;
  if     (fmode >=11 && fmode <= 19)neutrinoID= 14;
  else if(fmode >=21 && fmode <= 29)neutrinoID=-14;
  else if(fmode >=31 && fmode <= 39)neutrinoID= 12;
  else if(fmode >=41 && fmode <= 49)neutrinoID=-12;

  theVector->setDecayMode(fmode); 

  theVector->setNuType(neutrinoID);         
  

  theVector->setDecayParticleID(fppid); 

  theVector->setDecayMomentum(fppi*fnpi[0] * CLHEP::GeV, 
			      fppi*fnpi[1] * CLHEP::GeV, 
			      fppi*fnpi[2] * CLHEP::GeV);

  theVector->setDecayPos(fxpi[0] * CLHEP::cm, 
			 fxpi[1] * CLHEP::cm, 
			 fxpi[2] * CLHEP::cm);

  theVector->setCosToBeamAtDecay(fcospibm);

  theVector->setNuWeight(fnorm);
             

  theVector->setProductionMomentum(fppi0*fnpi0[0] * CLHEP::GeV, 
				   fppi0*fnpi0[1] * CLHEP::GeV, 
				   fppi0*fnpi0[2] * CLHEP::GeV);

  theVector->setProductionPos(fxpi0[0] * CLHEP::cm, 
			      fxpi0[1] * CLHEP::cm, 
			      fxpi0[2] * CLHEP::cm);

  theVector->setCosToBeamAtProduction(fcospi0bm);


  theVector->setNvtx0(fnvtx0);

  theVector->setNuRnu(frnu * CLHEP::cm);

  theVector->setIdfd(fidfd);

  // Translate if ND flux plane is not centered at ND280 origin (magnet center)
  theVector->setNuPos((fxnu + fFluxPlanePos[fidfd-1][0]) * CLHEP::cm, 
		      (fynu + fFluxPlanePos[fidfd-1][1]) * CLHEP::cm, 
		      (   0 + fFluxPlanePos[fidfd-1][2]) * CLHEP::cm); 

  theVector->setPrimaryParticleID(fgipart);
  theVector->setPrimaryParticlePos(fgpos0[0] * CLHEP::cm,
				   fgpos0[1] * CLHEP::cm, 
				   fgpos0[2] * CLHEP::cm);
  theVector->setPrimaryParticleMomentum(fgvec0[0]*fgamom0 * CLHEP::GeV,
					fgvec0[1]*fgamom0 * CLHEP::GeV,
					fgvec0[2]*fgamom0 * CLHEP::GeV);


  theVector->setNuVersion(version);

  // jnubeam > 2010d
  if (2010.35 < version) {
    theVector->setNumIntSteps(fng);

    theVector->setIntermediateParticleMomentum(fgpx, fgpy, fgpz);
    theVector->setCosToBeamOfIntermediate(fgcosbm);
    theVector->setIntermediateParticlePos(fgvx, fgvy, fgvz);
    theVector->setIntermediateParticleID(fgpid);
    theVector->setIntermediateParticleMec(fgmec);


    theVector->setNuEnergySK(fEnusk * CLHEP::GeV);
    theVector->setNuWeightSK(fnormsk);
    theVector->setNuWeightAccept(fanorm);
    
    //jnubeam >= 2011a
    if (2011 <= version) {
      theVector->setNuGmat(fgmat);
      theVector->setNuGdistc(fgdistc);
      theVector->setNuGdistal(fgdistal);
      theVector->setNuGdistti(fgdistti);
      theVector->setNuGdistfe(fgdistfe);
    }

    theVector->setNuTuneid(tuneid);
    theVector->setNuNtrig(fntrig);
    theVector->setNuPint(fpint);
    theVector->setNuBpos(fbpos) ;
    theVector->setNuBtilt(fbtilt) ;
    theVector->setNuBrms(fbrms) ;
    theVector->setNuEmit(femit) ;
    theVector->setNuAlpha(falpha) ;
    theVector->setNuHcur(fhcur) ;
    theVector->setNuRand(frand);
    theVector->setNuRseed(frseed);

  }

  return;
}

