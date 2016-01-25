// Event Rate calculation:
// For a given flux (jnubeam files) and ND280 detector geometry
// this program calculates the expected number of neutrino interactions
// as well as the maximum interaction probability for any given neutrino
// vector from jnubeam. It also produces a set of histograms which show
// the energy distribution for each species of neutrino, absolutely normalized
// for 10^{21} POT.
//
// The total number of expected interactions can be used to normalize
// generated neutrino samples to an equivalent POT
// The maximum probability is used to optimize the event generation process.
//    it is expressed in units of 10^{21} (i.e. the actual probabilty for the
//    neutrino with maximum interaction probability is the number reported
//     at the end divided by 10^{21}.

//
// Options
// -g geo.root specify geometry file (default is "nd280geometry.root")
//
// -v +V1-V2           include V1 in geometry, exclude -V2
//    
//
//  There are two ways to specify the jnubeam files (one of the two must be specified):
// -f jnubeamfile.root  run on a single file (called jnubeamfile.root)
// -s filestem n m      run on a series of files with the form filestem.i.root
//                      where i runs from n to m
//                      the default filestem is nu.nd280.i.root
//
// -p NEUTRINO          Select neutrino species to run 0 = all is default
//                      +/-12 for nue(bar), +/- 14 for numu(bar)
//
// -t treename          This specifies the tree containing the neutrino vectors 
//                      h3002 is default
//
// -d detector          This specifies "ND setting", i.e. which patch of off-axis
//                      solid angle to use. The default for ND280 is 5 (Basket). 
//                      For INGRID 3/4, use 3 only. 13 for sand simulation.
//
// -o out.root output   root file name (default will have the form 
//                      filestem.n-m_neutrate.root (in the original flux file dir)


#include "SystemOfUnits.h"
#include "TObjString.h"
#include "TNuFlux.h"
#include "TNuVector.h"
#include "TNuVertex.h"
#include "TNuTrajectory.h"
#include "TVector3.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGraph.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1D.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h" 
#include "Exclude.h"
#include "TVectorD.h"
#include "TKey.h"
#include <iostream>  

int getArgs(int argc, char* argv[]);
TString geomFile     = "nd280geometry.root";
std::string  genVol       = "+Basket";
TString jnuFile      = "";
TString jnuDir       = "flux";
TString jnuStem      = "nu.nd280.";
TString jnuTreeName  = "h3002";
TString jnuVerTreeName = "h1000";
TString jnuNDTreeName  = "h3000";
Int_t   jnuStart     = -1;
Int_t   jnuEnd       = -1;
Int_t   jnuDetector  = 5;
Int_t   nuTypePDG    = 0;
Int_t   nuTypeJnuBeam = 0;
Bool_t  singleFile   = true;
TString outputFile   = "";
std::string nuGenerator   = "gibuu"; //"neut";  // TF: NEW!!!

int main(int argc, char *argv[]){

  // process the arguments
  int args = getArgs(argc, argv);
  if(args != 0){
    std::cerr << "Usage " << std::endl;
    return 0;
  }

  // set the output file name if blank
  if(outputFile==""){
    // if jnuFile is blank, we are using the standard jnubeam production
    if(jnuFile==""){
      outputFile = Form("%s.%d-%d_neutrate.root", 
			jnuStem.Data(),	jnuStart, jnuEnd);
    } else {

      std::string filename = jnuFile.Data();
      /*
      // Get basename of flux file
      std::string::size_type start_pos = filename.rfind("/");
      if (start_pos == std::string::npos) start_pos = 0; else ++start_pos;
      std::string basename = filename.substr(start_pos);
      /**/

      std::string::size_type ext_pos = filename.rfind(".root");
      if (ext_pos == std::string::npos) ext_pos = 0; else ++ext_pos;
      filename = filename.substr(0,ext_pos-1);

      outputFile = Form("%s_neutrate.root", filename.c_str());

      std::cout << "Output file name: " << outputFile << std::endl;
    }
  }
  

  // open the geometry file
  TFile fnd280geom(geomFile.Data());
  
  // Get TGeoManager key name
  TIter next(fnd280geom.GetListOfKeys());
  TKey *key;
  Int_t foundGeo=0;
  TString geoKeyName, className;
  while ((key=(TKey*)next())) {
    geoKeyName   = key->GetName();
    className = key->GetClassName();

    if (!className.CompareTo("TGeoManager") && geoKeyName.BeginsWith("ND280Geometry")) {
      foundGeo=1;
      break;
    }
  }
  if (!foundGeo) {
    std::cerr << "Error: Could not find TGeoManager object in " << geomFile << std::endl;
    exit (-1);
  }
  
  std::cout << "Geometry file " << geomFile  << std::endl;
  std::cout << "Geometry hash " << geoKeyName  << std::endl;
  std::cout << "Volumes       " << genVol    << std::endl;
  if (singleFile)
    std::cout << "jnubeam file  " << jnuFile   << std::endl;
  else
    std::cout << "jnubeam stem start/end  "    << jnuStem << " " 
	      << jnuStart << " " << jnuEnd << std::endl;



  // Get the geometry manager
  TGeoManager *geom = (TGeoManager*)fnd280geom.Get(geoKeyName);


  TFile fout(outputFile.Data(), "RECREATE");

  TH1D* hnu       = new TH1D("hnu"      , ";E_{#nu} (GeV); Events/0.10 GeV", 400, 0.0, 20.0);
  TH1D* hev       = new TH1D("hev"      , ";E_{#nu} (GeV); Events/0.10 GeV", 400, 0.0, 20.0);
  TH1D* hev_numu  = new TH1D("hev_numu" , ";E_{#nu} (GeV); Events/0.10 GeV", 400, 0.0, 20.0);
  TH1D* hev_numub = new TH1D("hev_numub", ";E_{#nu} (GeV); Events/0.10 GeV", 400, 0.0, 20.0);
  TH1D* hev_nue   = new TH1D("hev_nue"  , ";E_{#nu} (GeV); Events/0.10 GeV", 400, 0.0, 20.0);
  TH1D* hev_nueb  = new TH1D("hev_nueb" , ";E_{#nu} (GeV); Events/0.10 GeV", 400, 0.0, 20.0);
  TH1D* hev_prob  = new TH1D("hev_prob" , ";Probability; Events",        1000, 0.0, 10);
  TH1D* hev_int  = new TH1D("hev_int" , ";Probability; Events/10",        5000, 0.0, 50000);


  TTree* theProbTree = new TTree("pTree", "Probabilities");
  TObjString* fNuFileName = 0;
  Int_t ev;
  Double_t prob, nuprob, beamwgt;
 
  theProbTree->Branch("NuFileName",    "TObjString", &fNuFileName,       32000, 1   );
  theProbTree->Branch("ev"      ,&ev     , "ev/I"          );
  theProbTree->Branch("prob"    ,&prob   , "prob/D"        );
  theProbTree->Branch("nuprob"  ,&nuprob , "nuprob/D"      );
  theProbTree->Branch("beamwgt" ,&beamwgt, "beamwgt/D"     );

  TGeoVolume * theTopVol = geom->GetTopVolume();
  PrepareVolume(genVol, theTopVol);  
  
  TNuTrajectory theNuTraj;
  
  theNuTraj.SetGeoManager(geom);
  theNuTraj.SetRandomSeed(0);
  theNuTraj.SetGlobalVol(theTopVol);
  theNuTraj.SetXsecFunc(nuGenerator);  // TF: NEW and works for all generators


  // Sand flux starts way back
  if (jnuDetector==13)
    theNuTraj.SetStartZ(-5000 * CLHEP::cm);
  else
    theNuTraj.SetStartZ(-1000 * CLHEP::cm);

  Int_t nFluxFiles = (abs(jnuEnd - jnuStart) + 1);
  TNuFlux* theFluxFiles;
  if(singleFile){
    theFluxFiles = new TNuFlux(jnuFile,jnuTreeName,jnuVerTreeName,jnuNDTreeName);
    theFluxFiles->addToChain(jnuFile);
    if (theFluxFiles->checkVerTree(jnuFile))
      return -1;
  } else {
    theFluxFiles = new TNuFlux(Form("%s.%d.root", jnuStem.Data(), jnuStart),jnuTreeName.Data(),jnuVerTreeName.Data(),jnuNDTreeName.Data());
    for(int i =  jnuStart; i <= jnuEnd; i++) {
      theFluxFiles->addToChain(Form("%s.%d.root",  jnuStem.Data(), i));
      if (theFluxFiles->checkVerTree(Form("%s.%d.root", jnuStem.Data(), i)))
	return -1;
    }
  }
  Int_t nChainFluxFiles = theFluxFiles->getNChainFluxFiles();
  Int_t nHaddFluxFiles = theFluxFiles->getNHaddFluxFiles();

  if (nFluxFiles != nChainFluxFiles) {
    std::cout << "Error in chaining flux files: " << nFluxFiles << " != " << nChainFluxFiles << std::endl;
    return -1;
  }    

  theFluxFiles->setDetector(jnuDetector);
  theFluxFiles->setNuType(nuTypeJnuBeam);
  std::cout << "Total events " << theFluxFiles->getNEntries() << " amongst " << nHaddFluxFiles << " hadd'ed flux files of 10^21 POT each, from " << nFluxFiles << " chained input files." << std::endl;
  

  //                        all  numu  numubar nue nuebar 
  Double_t maxIntRate[5] = {0,   0,    0,      0,      0    };
  Double_t nExpInts[5]   = {0,   0,    0,      0,      0    };
  TNuVector* theNuVector = 0;
  Int_t iEvent = 0; 

  while(theNuVector = theFluxFiles->getVector()){

    if(iEvent%1000 == 0)std::cout << "Vector " << iEvent << std::endl;

    theNuTraj.Swim(theNuVector);
    theNuTraj.ExtractNuclei();
    theNuTraj.CalcXsecs();
    theNuTraj.CalcScatCoefficients();
    theNuTraj.CalcSumIntProb();

    Double_t sumProb    = theNuTraj.GetSumIntProb();
    Double_t beamWeight = theNuVector->getNuWeight(); 
    Double_t energy     = theNuVector->getNuEnergy();
    Int_t    nuType     = theNuVector->getNuType();

    Double_t intPer1021    = sumProb*beamWeight;

    //cout << sumProb << " " << beamWeight << " " << intPer1021 << " " << energy << " " << nuType << endl;

    // fill ntuple variables
    if(fNuFileName) { 
      delete fNuFileName;
      fNuFileName = 0;
    }
    fNuFileName = new TObjString((theNuVector->getNuFileName()).Data());
    ev      = theNuVector->getNuFluxEntry();
    prob    = intPer1021;
    nuprob  = sumProb;
    beamwgt = beamWeight;
    theProbTree->Fill();


    if( intPer1021 >= 0 && 
        ( (jnuDetector!=13 && intPer1021 < 1e10) || 
	  (jnuDetector==13 && intPer1021 < 1e15) ) ) 
      {

      // Number of interactions 

      // all neutrinos
      nExpInts[0]     += intPer1021;

      // for each species
      if(       nuType ==  14 ){ // numu
	nExpInts[1] += intPer1021;
      } else if(nuType == -14 ){ // numubar
	nExpInts[2] += intPer1021; 
      } else if(nuType ==  12 ){ // nue
	nExpInts[3] += intPer1021;
      } else if(nuType == -12 ){ // nuebar
	nExpInts[4] += intPer1021;
      }

      // Maximum interaction probabilities

      // all neutrinos
      if (intPer1021 > maxIntRate[0]) maxIntRate[0] = intPer1021;

      // for each species
      if(       nuType ==  14 && intPer1021 > maxIntRate[1]){ // numu
	maxIntRate[1] = intPer1021;
      } else if(nuType == -14 && intPer1021 > maxIntRate[2]){ // numubar
	maxIntRate[2] = intPer1021; 
      } else if(nuType ==  12 && intPer1021 > maxIntRate[3]){ // nue
	maxIntRate[3] = intPer1021;
      } else if(nuType == -12 && intPer1021 > maxIntRate[4]){ // nuebar
	maxIntRate[4] = intPer1021;
      }

    } else {
      std::cerr << "Negative or large interaction rate reported: " << intPer1021 <<std::endl;
      std::cerr << "in flux file: " << (theNuVector->getNuFileName()).Data() << ", vector: " << ev << std::endl;
      std::cerr << "Ignoring in determing total Interaction rate" << std::endl;
    }
    if(sumProb>0)hev_prob->Fill(sumProb*1e13);
    hev_int->Fill(intPer1021);

    hnu->Fill(energy/CLHEP::GeV);
    hev->Fill(energy/CLHEP::GeV, intPer1021);
    if     (nuType == 14)hev_numu ->Fill(energy/CLHEP::GeV, intPer1021);
    else if(nuType ==-14)hev_numub->Fill(energy/CLHEP::GeV, intPer1021);
    else if(nuType == 12)hev_nue  ->Fill(energy/CLHEP::GeV, intPer1021);
    else if(nuType ==-12)hev_nueb ->Fill(energy/CLHEP::GeV, intPer1021);

    delete theNuVector;
    iEvent++;

    

  }

  theNuTraj.Clear();

  for (int inu=0; inu<5; inu++)
    nExpInts[inu] = nExpInts[inu]/(double)nHaddFluxFiles;
  

  std::cout << std::endl << "NuType  maxIntProb   nExpInts/10^21 POT" << std::endl;
  std::cout << "all     " << maxIntRate[0] <<  "   " << nExpInts[0] << std::endl;
  std::cout << "numu    " << maxIntRate[1] <<  "   " << nExpInts[1] << std::endl;
  std::cout << "numubar " << maxIntRate[2] <<  "   " << nExpInts[2] << std::endl;
  std::cout << "nue     " << maxIntRate[3] <<  "   " << nExpInts[3] << std::endl;
  std::cout << "nuebar  " << maxIntRate[4] <<  "   " << nExpInts[4] << std::endl;

  std::cout << std::endl << "Stored in file: " << outputFile.Data() << std::endl << std::endl;

  TVectorT<double> theMaxIntRate(5, maxIntRate);
  TVectorT<double> theNexpInts(5, nExpInts);
  fout.cd();
  theMaxIntRate.Write("MaxIntRate");
  theNexpInts.Write("nExpInts");
  fout.Write(); fout.Close();
  
  return 0;
  
}



int getArgs(int argc, char* argv[]){

  while( (argc > 1) && (argv[1][0] == '-') ){
    switch(argv[1][1]){

    case 'd':
      jnuDetector    = atoi(argv[2]);
      ++argv; --argc;
      break;
    case 'f':
      if(singleFile){
	jnuFile        = argv[2];
	++argv; --argc;
	break;
      } else {
	std::cerr << "Don't specify a jnubeam file if using the stem/range option" << std::endl;
	return -1;
      }
    case 'g': 
      geomFile = argv[2];
      ++argv; --argc;
      break;
    case 'o':
      outputFile     = argv[2];
      ++argv; --argc;
      break;
    case 'p':
      nuTypePDG  = atoi(argv[2]);
      if(nuTypePDG == 0){
	std::cout << "Running on all neutrino species " << std::endl;
	nuTypeJnuBeam = 0;
	/* Deprecated (issues with index matching in genev.cc)
      } else if (nuTypePDG == 12){
	std::cout << "Running nue interactions only " << std::endl;
	nuTypeJnuBeam = 3;
      } else if (nuTypePDG == -12){
	std::cout << "Running nuebar interactions only " << std::endl;
	nuTypeJnuBeam = 4;
      } else if (nuTypePDG == 14){
	std::cout << "Running numu interactions only " << std::endl;
	nuTypeJnuBeam = 1;
      } else if (nuTypePDG == -14){
	std::cout << "Running numubar interactions only " << std::endl;
	nuTypeJnuBeam = 2;
	*/
      } else {
	std::cout << "Warning: -p option is deprecated." << std::endl;
	std::cout << "Setting nuType = 0 and generating all species " << std::endl;
	nuTypeJnuBeam = 0;
      }
      ++argv; --argc;
      break;
    case 's':
      jnuStem = argv[2];
      ++argv; --argc;
      if(argv[2][0] != '-'){
	jnuStart    = atoi(argv[2]);
	++argv; --argc;
      }	else {
	std::cerr << "Need to specify starting index of jnubeam files " << std::endl;  
	return -1;
      }
      if(argv[2][0] != '-'){
	jnuEnd    = atoi(argv[2]);
	++argv; --argc;
      }	else {
	std::cerr << "Need to specify ending index of jnubeam files " << std::endl;  
	return -1;
      }
      singleFile = false;
      break;
    case 't':
      jnuTreeName    = argv[2];
      ++argv; --argc;
      break;
    case 'u':
      jnuVerTreeName    = argv[2];
      ++argv; --argc;
      break;
    case 'v':
      genVol  = argv[2];
      ++argv; --argc; 
      break;
    // TF: NEW OPTION: Generator
    case 'z':
      nuGenerator = argv[2];
      ++argv; --argc; 
      break;
    }
    
    ++argv; --argc;
  }

  return 0;

}
