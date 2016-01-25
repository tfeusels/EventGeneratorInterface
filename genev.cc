#include "SystemOfUnits.h"
#include "TNuFlux.h"
#include "TNuVector.h"
#include "TNuVertex.h"
#include "TNuTrajectory.h"


#include "TNeutOutput.h"


#include "TVector3.h"
#include "TLorentzVector.h"
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
#include "TChain.h"
#include <iostream>
//#include <algorithm>
#include "assert.h"
#include "necardC.h"


//
// Options
// -g geo.root specify geometry file (default is "nd280geometry.root")

// -v +V1-V2           include V1 in geometry, exclude -V2
//    

//  There are two ways to specify the jnubeam files 
//    (one of the two must be specified):
// -j jnubeamfile.root  run on a single file (called jnubeamfile.root)
// -s filestem n m      run on a series of files with the form filestem.[i].root
//                      where i runs from n to m
//                      the default filestem is nu.nd280.[i].root.
//                      This requires event_rate to be run on all corresponding
//                      flux files with default file name: filestem.[i]_neutrate.root

//  Other options
// Number of events/pot: one should specify only one of these.
// -n number of events  number of events to generate. 
//                      0 = exhaust available jnubeam vectors is default
//
// -e effective number of protons on target to generate (in units of 10^{21})                         
// -N number of expected events/10^21 POT.
//
//
// -r random number seed Default = 0 (TRandom3 setting).
// -c RECYCLE           rethrow interaction RECYCLE number of times.
//                      default RECYCLE = 1
// -w                   rewind the jnubeam tree if the vectors are exhausted and the
//                      specified number of events/pot has not been generated. Default is false
//                       i.e. no -w option
//
// -p NEUTRINO          Select neutrino species to run 0 = all is default
//                      +/-12 for nue(bar), +/- 14 for numu(bar)
//
// -t treename          This specifies the tree containing the neutrino 
//                      vectors. h3001 is default
// -d detector          This specifies "ND setting", i.e. which patch of 
//                      solid angle to use. 
//                      The default for ND280 is 5 (Basket).
//                      6: ND280 Magnet
//                      For INGRID 3/4, use 3 only since ntuple contains 
//                      both 3 and 4. 
//                      13: Sand/Dirt/Muck
//
// -o out.root output   root file name (default will have the form 
//                      eventrate_filestem_n-m_+VOLUME-VOLUME.root
//
// -f neut/rootracker   format of output, neut (Deprecated) or rootracker
//
// -m INTRATE           maximum interaction probability for vector
//                      this is used to optimize the rejeciton method
//                      for determining whether a given neutrino interacts.
//                      It is calcualted for a given geometry/flux files
//                      by event_rate and is by default 0 and you must set
//                      it by this option or with input event_rate file
//
// -i int.pro file      file containg precalculated interaction probabilities for
//                      each neutrino vector. If this is not supplied, the probabilities
//                      are calculated on the fly
//
// -q random order      randomize order of chained flux files (default = 0: off)
// -Q random start      random starting entry in chain or input file (default = 0: off)


int getArgs(int argc, char* argv[]);
TString  geomFile          = "nd280geometry.root";
string   genVol            = "+Basket";
TString  jnuFile           = "";
TString  jnuDir            = "flux";
TString  jnuStem           = "nu.nd280.";
TString  jnuTreeName       = "h3002";
TString  jnuVerTreeName    = "h1000";
TString  jnuNDTreeName     = "h3000";
Int_t    jnuStart          =  -1;
Int_t    jnuEnd            =  -1;
Int_t    jnuDetector       =   5;
Int_t    nuTypePDG         =   0;
Int_t    nuTypeJnuBeam     =   0;
Int_t    nRecycle          =   1;
Bool_t   rewindJnuBeam     = false;
Int_t    randomSeed        =   0;
Int_t    nEventsToGenerate =   0;
Double_t potToGenerate     =   0;
Double_t eventsPer1021     =   0;
Bool_t   singleFile        = true;
Int_t    numberOfEvents    =   0;
Bool_t   repeatFiles       = false;
TString  outputFile        = "";
TString  outputFileFormat  = "rootracker";
Double_t maxIntProb        = 0;
TString  intProbFile       = "";
Bool_t   rndJnuOrder       = false;
Bool_t   rndEntryStart     = false;
std::string nuGenerator    = "neut";

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
      outputFile = Form("%s.%d-%d_neutgenev.root",jnuStem.Data(), jnuStart, jnuEnd);
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

      outputFile = Form("%s_neutgenev.root", filename.c_str());

      std::cout << "Output file name: " << outputFile << std::endl;
    }
  }

  
  TRandom3 theRandom; theRandom.SetSeed(randomSeed);
  

  // Shuffle flux file indices
  const int nFluxFiles = abs(jnuEnd-jnuStart)+1;
  int fluxFileIndex[nFluxFiles];
  
  if (!singleFile) {
    for (int i=jnuStart; i<=jnuEnd; i++) {
      fluxFileIndex[i-jnuStart] = i;  // fill the array in order
    }

    // Randomize order
    if (rndJnuOrder) {
      for (int i=0; i<(nFluxFiles-1); i++) {
	int r = i + (theRandom.Rndm() * (nFluxFiles-i)); // Random remaining position.
	int temp = fluxFileIndex[i]; 
	fluxFileIndex[i] = fluxFileIndex[r]; 
	fluxFileIndex[r] = temp;
      }
      
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

  // get the geometry manager
  TGeoManager *geom = (TGeoManager*)fnd280geom.Get(geoKeyName);

  TGeoVolume * theTopVol = geom->GetTopVolume();
  PrepareVolume(genVol, theTopVol);  

  // Output information to screen
  std::cout << std::endl << "Geometry file " << geomFile  << std::endl;
  std::cout << "Geometry hash " << geoKeyName  << std::endl;
  std::cout << "Volumes       " << genVol    << std::endl << std::endl;



  // get the interaction probability tree if specified
  Bool_t calcProbOnTheFly = true; // by default calculate the probabilites on the fly
  TChain *theIntProbTree = new TChain("pTree");
  TObjString* NuFileName = NULL; 
  Int_t ev; 
  Double_t prob, nuprob, nuwgt;
  TBranch *b_NuFileName, *b_ev, *b_prob;
  TVectorT<double> *theMaxIntRate = NULL;
  TVectorT<double> *theNexpInts = NULL;
  double *maxIntRate;
  double *nExpInts;

  if( intProbFile != "" ){


    std::cout << "Loading max interaction rate, expected # of" << std::endl;
    std::cout << "interactions and probability trees from file(s):" << std::endl;
    
    TString theIntProbFileName = intProbFile.Data(); 


    // Warning: Following lines are commented out assuming input setup files are NOT chained
    //          i.e. produced with event_rate with a hadd'ed (or TChain'ed) flux file(s)
    //for (int ifile=0; ifile<nFluxFiles; ifile++) {
      //if (!singleFile) theIntProbFileName = Form("%s.%d_neutrate.root",jnuStem.Data(),fluxFileIndex[ifile]); // Remember to change below too if modifying
      
      std::cout << theIntProbFileName.Data() << std::endl;
      
      // Load pre-determined maximum interaction rate and expected event rates
      TFile* fint = new TFile(theIntProbFileName.Data(), "READ");

      theMaxIntRate = (TVectorT<double>*)fint->Get("MaxIntRate");
      theNexpInts = (TVectorT<double>*)fint->Get("nExpInts");
      if (theMaxIntRate == NULL){
	std::cerr << "Error: Maximum Interaction Probability Vector not found " << std::endl;
	exit (-1);
      }
      else if (theNexpInts == NULL){
	std::cerr << "Error: Number of Expected Events Vector not found " << std::endl;
	exit (-1);
      }
      else {
	maxIntRate = theMaxIntRate->GetMatrixArray();
	nExpInts = theNexpInts->GetMatrixArray();

	// Use maximum interaction probability for specific neutrino type
	if(nuTypePDG == 0){
	  maxIntProb = std::max(maxIntProb,maxIntRate[0]);
	  eventsPer1021 += nExpInts[0];

	} else if (nuTypePDG == 14){
	  maxIntProb = std::max(maxIntProb,maxIntRate[1]);
	  eventsPer1021 += nExpInts[1];

	} else if (nuTypePDG == -14){
	  maxIntProb = std::max(maxIntProb,maxIntRate[2]);
	  eventsPer1021 += nExpInts[2];

	} else if (nuTypePDG == 12){
	  maxIntProb = std::max(maxIntProb,maxIntRate[3]);
	  eventsPer1021 += nExpInts[3];

	} else if (nuTypePDG == -12){
	  maxIntProb = std::max(maxIntProb,maxIntRate[4]);
	  eventsPer1021 += nExpInts[4];

	} else { // Default to all neutrinos
	  maxIntProb = std::max(maxIntProb,maxIntRate[0]);
	  eventsPer1021 += nExpInts[0];
	}
      }
      fint->Close();
      //}
      //eventsPer1021 = (Double_t)eventsPer1021/(Double_t)nFluxFiles;

    std::cout << std::endl << "Using Max Interaction Probability from file(s) = " << maxIntProb << std::endl;
    std::cout << "Using Expected Event Rate from file(s) = " << eventsPer1021 << " / 10^21 POT" << std::endl << std::endl;


    // Loaded precaculated probability tree
    theIntProbTree->SetBranchAddress("NuFileName", &NuFileName, &b_NuFileName);       
    theIntProbTree->SetBranchAddress("ev",     &ev, &b_ev);
    theIntProbTree->SetBranchAddress("prob",   &prob, &b_prob);

    // Warning: Following lines are commented out assuming input setup files are NOT chained
    //          i.e. produced with event_rate with a hadd'ed (or TChain'ed) flux file(s)
    //for (int ifile=0; ifile<nFluxFiles; ifile++) {

      //if (!singleFile) theIntProbFileName = Form("%s.%d_neutrate.root",jnuStem.Data(),fluxFileIndex[ifile]); // Remember to change above too if modifying    

      if (!theIntProbTree->Add(Form("%s/pTree",theIntProbFileName.Data()))) {
	std::cerr << "File or tree not found at: " << theIntProbFileName.Data() << std::endl;
	exit (-1);
      }
      //}
    
    calcProbOnTheFly = false;
    
  } else {
    std::cerr << "You need to specify a setup file generated with event_rate via command line argument: -i " << std::endl;
    exit (-1);
  }
 





  // Prepare the trajectory/rate calculations
  TNuTrajectory theNuTraj;
  theNuTraj.SetGeoManager(geom);
  theNuTraj.SetRandomSeed(randomSeed);
  theNuTraj.SetGlobalVol(theTopVol);
  theNuTraj.SetXsecFunc(nuGenerator); //TF: NEW! WORKS FOR ANY GENERATOR

  // Sand flux starts way back
  if (jnuDetector==13)
    theNuTraj.SetStartZ(-5000 * CLHEP::cm);
  else
    theNuTraj.SetStartZ(-1000 * CLHEP::cm);

  // The Neutrino event generator
  //TNeutEventGenerator* theVectorGenerator = new TNeutEventGenerator();
  TNuEventGenerator* theVectorGenerator;
  if(nuGenerator == "neut")
    theVectorGenerator = new TNeutEventGenerator();


  // Prepare the input flux files
  TNuFlux* theFluxFiles;
  

  if (singleFile)
    std::cout << "jnubeam single file  " << jnuFile   << std::endl;
  else {
    std::cout << "jnubeam chain stem start/end  "    << jnuStem << " " 
            << jnuStart << " " << jnuEnd << std::endl;
    if (rndJnuOrder)
      std::cout << "Flux files chained in random order." << std::endl;
    else
      std::cout << "Flux files chained sequentially." << std::endl;
  }

  // Single file
  if(singleFile){
    theFluxFiles = new TNuFlux(jnuFile,jnuTreeName,jnuVerTreeName,jnuNDTreeName);
    theFluxFiles->addToChain(jnuFile);
    if (theFluxFiles->checkVerTree(jnuFile)) {
      std::cerr << "GENEV aborting (single file)" << std::endl;
      return -1;
    }

  // Chain of files
  } else {
    
    if (rndJnuOrder) {
      theFluxFiles = new TNuFlux(Form("%s.%d.root", jnuStem.Data(), fluxFileIndex[0]),jnuTreeName.Data(),jnuVerTreeName.Data(),jnuNDTreeName.Data());
      for(int i = 0; i < nFluxFiles; i++) {
	TString fluxFileName = Form("%s.%d.root", jnuStem.Data(), fluxFileIndex[i]);
	theFluxFiles->addToChain(fluxFileName.Data());
	if (theFluxFiles->checkVerTree(fluxFileName.Data())) {
	  std::cerr << "GENEV aborting (chained files)" << std::endl;
	  return -1;
	}
	std::cout << "Flux file chained: " << fluxFileName.Data() << std::endl;
      }

    } else { // Sequential chaining
    
      theFluxFiles = new TNuFlux(Form("%s.%d.root", jnuStem.Data(), jnuStart),jnuTreeName.Data(),jnuVerTreeName.Data(),jnuNDTreeName.Data());
      for(int i = jnuStart; i <= jnuEnd; i++) {
	theFluxFiles->addToChain(Form("%s.%d.root", jnuStem.Data(), i));
	if (theFluxFiles->checkVerTree(Form("%s.%d.root", jnuStem.Data(), i))) {
	  std::cerr << "GENEV aborting" << std::endl;
	  return -1;
	}
      }
    }
  }
  Int_t nChainFluxFiles = theFluxFiles->getNChainFluxFiles();
  Int_t nHaddFluxFiles = theFluxFiles->getNHaddFluxFiles();

  if (nFluxFiles != nChainFluxFiles) {
    std::cout << "Error in chaining flux files: " << nFluxFiles << " != " << nChainFluxFiles << std::endl;
    return -1;
  }
  if (!singleFile && nHaddFluxFiles>nChainFluxFiles) {
    std::cout << "Error: Do not chain hadd'ed flux files!" << std::endl;
    return -1;
  }

  theFluxFiles->setDetector(jnuDetector);
  theFluxFiles->setNuType(nuTypeJnuBeam);
  theFluxFiles->setRewind(rewindJnuBeam);

  std::cout << "Total JnuBeam Vectors " << theFluxFiles->getNEntries() << std::endl ;

  
  // Get random starting position
  Int_t startingEntry=0;
  Int_t currSkip=0;;
  if (rndEntryStart) {
    startingEntry = TMath::Nint( theRandom.Rndm() * (theFluxFiles->getNEntries()-1) );
    theFluxFiles->setVectorIndex(startingEntry);
    std::cout << "Skipping the first " << startingEntry << " vectors." << std::endl << std::endl;
  }
  
  // See if we want to generate according to POT;
  if(potToGenerate > 0){
    if(eventsPer1021 >0){
      Double_t meanNEvents = eventsPer1021 * potToGenerate;
      nEventsToGenerate = theRandom.Poisson(meanNEvents);
    } else {
      std::cerr << "To generate POT equivalent samples, both the POT to "
		<< "generate and the number of events/10^21 POT must be specified"
		<< std::endl;
      exit (-1);
    }
  }

  // the output file
  TNeutOutput* theNeutOutput;
  if       (outputFileFormat == "neut"      ){
    //theNeutOutput = new TNeutOutput(outputFile.Data(), 0, theFluxFiles->getVersion());
    std::cerr << "NEUT output format is no longer supported. Use rootracker format instead." << std::endl;
    return -1;
  } else if(outputFileFormat == "rootracker"){
    theNeutOutput = new TNeutOutput(outputFile.Data(), 1, theFluxFiles->getVersion());
  } else {
    std::cerr << "Unknown format " << outputFileFormat << std::endl;
    return -1;
  }

  theNeutOutput->SetTreeWeight(potToGenerate*1e21);

  std::cout << std::endl << "Output file   " << outputFile << std::endl << std::endl;


  TNuVector* theNuVector       =     0;
  Int_t validJnuBeam           =     0;  
  Int_t generatedInteractions  =     0;  

  if (!maxIntProb) {
    std::cerr << "You must set the maximum interaction rate via -m or input event_rate file(s)" << std::endl;
    exit (-1);
  }
  theNuTraj.SetMaxIntProb(maxIntProb);


  // The event loop
  while( (theNuVector = theFluxFiles->getVector()) &&
	 (nEventsToGenerate == 0 || generatedInteractions < nEventsToGenerate) ){

    //if(validJnuBeam%1000 == 0)std::cout << "JnuBeam Vector " << validJnuBeam << std::endl;

    // if we have precalculated the interaction probability tree, get it here
    Double_t probability = 0;
    int event;
    TString fluxfilename;
    int NuFluxEntry;
    
    if(!calcProbOnTheFly){

      if(NuFileName) { 
	delete NuFileName;
	NuFileName = 0;
      }

      event = theNuVector->getEntry();
      theIntProbTree->GetEntry(event);

      TString fluxfilename_prob = NuFileName->GetString();
      
      fluxfilename = theNuVector->getNuFileName();      
      NuFluxEntry = theNuVector->getNuFluxEntry();
      
      //cout << event << " | " << fluxfilename_prob.Data() << ":" << fluxfilename.Data() << " | " <<  ev << ":" << NuFluxEntry <<  endl;
      
      // Make sure flux file entry matches probability entry
      assert(!fluxfilename.CompareTo(fluxfilename_prob));
      assert(NuFluxEntry == ev);

      //cout << fluxfilename.Data() << " " << fluxfilename_prob.Data() << " " << NuFluxEntry << " " << ev << endl;

      probability = prob/maxIntProb;
    }

    if(calcProbOnTheFly || theRandom.Rndm() < probability){
      theNuTraj.Swim(theNuVector);
      theNuTraj.ExtractNuclei();
      theNuTraj.CalcXsecs();
      theNuTraj.CalcScatCoefficients();
      theNuTraj.CalcSumIntProb();
      theNuTraj.SetBeamWeight(theNuVector->getNuWeight());

      for(int idraw = 0; idraw < nRecycle; idraw++){

	if(!calcProbOnTheFly || theNuTraj.DrawInteraction()){
	
	  Double_t energy = theNuVector->getNuEnergy();
	  
	  TNuVertex* theNuVertex = theNuTraj.DrawVertexFromPath();
	  TLorentzVector theVertex = theNuVertex->getVertex();

	  // generate an event;
	  float* pnu = new float[3];
	  TVector3 theNuDir = theNuVector->getNuDir();
	  pnu[0] = (float)(energy/CLHEP::GeV * theNuDir.X());
	  pnu[1] = (float)(energy/CLHEP::GeV * theNuDir.Y());
	  pnu[2] = (float)(energy/CLHEP::GeV * theNuDir.Z());
	  
	  Int_t nuType = theNuVector->getNuType();
	  
	  // Error code output from neutev do determine if we should skip 
	  // to next flux vector
	  Int_t ierr = 1;
	  
	  Int_t nTries = 0;
	  while(ierr == 1) {
	    theVectorGenerator->getVector(theNuVector->getNuType(),
					  pnu, 
					  theNuVertex->getZ(), 
					  theNuVertex->getIA(),
					  ierr);
	    nTries++;
	    if(nTries > 10){
	      break;
	    }
	  }
	  delete pnu;

	  // Break out of draw loop (and skip this flux vector)
	  if (ierr == 1) {
	    std::cerr << "genev: No interaction in neutev after 10 tries. skipping " << event << " | " << " Flux file:" << fluxfilename.Data() << " | " <<  ev << ":" << NuFluxEntry << std::endl;
	    delete theNuVertex;
	    break;

          // Good interaction
	  } else {
	    theNeutOutput->FillTree(generatedInteractions, theNuVector, 
				    theNuVertex, theVectorGenerator);
	    generatedInteractions++;
	    
	    if (neutcard_.quiet==0)
	      if(generatedInteractions%1000 == 0)std::cout << "Number of generated interactions " << generatedInteractions << std::endl;
	  }

	  delete theNuVertex;
	  
	} // generate the interaction
      } // loop over draws
    } // draw based on precalculated probability or go ahead with on the fly calculation

    delete theNuVector;
    validJnuBeam++;
    
  }

  std::cout << "Generated " << generatedInteractions 
	    << " in " << validJnuBeam << " valid Jnubeam vectors. " 
	    << std::endl;
  theNuTraj.Clear();
  theNeutOutput->WriteTree();
  theNeutOutput->CloseFile();

  return 0;
  
}




int getArgs(int argc, char* argv[]){

  while( (argc > 1) && (argv[1][0] == '-') ){

    switch(argv[1][1]){
    case 'c':
      nRecycle    = atoi(argv[2]);
      ++argv; --argc;
      break;
    case 'd':
      jnuDetector    = atoi(argv[2]);
      ++argv; --argc;
      break;
    case 'f':
      outputFileFormat = argv[2];
      ++argv; --argc;
      break;
    case 'g': 
      geomFile = argv[2];
      ++argv; --argc;
      break;
    case 'j':
      if(singleFile){
        jnuFile        = argv[2];
        ++argv; --argc;
        break;
      } else {
	std::cerr << "Don't specify a jnubeam file if using the stem/range option" << std::endl;
        return -1;
      }
    case 'm':
      maxIntProb = atof(argv[2]);
      ++argv; --argc;
      break;
    case 'n':
      nEventsToGenerate = atoi(argv[2]);
      ++argv; --argc;
      if(potToGenerate > 0){
	std::cerr << "You've specified both number of events to generate and "
		  << "POT to generate: will use number of events " << std::endl;
      }
      break;
    case 'e':
      potToGenerate = atof(argv[2]);
      ++argv; --argc;
      if(nEventsToGenerate > 0){
	std::cerr << "You've specified both number of events to generate and "
		  << "POT to generate: will use number of events " << std::endl;
      }
      break;
    case 'i':
      intProbFile  = argv[2];
      ++argv; --argc;
      break;
    case 'N':
      eventsPer1021 = atof(argv[2]);
      ++argv; --argc;
      break;
    case 'w':
      rewindJnuBeam = atoi(argv[2]);
      ++argv; --argc;
      break;
    case 'o':
      outputFile     = argv[2];
      //std::cout << outputFile << std::endl;
      ++argv; --argc;
      break;
    case 'p':
      nuTypePDG  = atoi(argv[2]);
      if(nuTypePDG == 0){
	std::cout << "Running on all neutrino species " << std::endl;
	nuTypeJnuBeam = 0;
      } else if (nuTypePDG == 14){
	std::cout << "Running numu interactions only " << std::endl;
	nuTypeJnuBeam = 1;
      } else if (nuTypePDG == -14){
	std::cout << "Running numubar interactions only " << std::endl;
	nuTypeJnuBeam = 2;
      } else if (nuTypePDG == 12){
	std::cout << "Running nue interactions only " << std::endl;
	nuTypeJnuBeam = 3;
      } else if (nuTypePDG == -12){
	std::cout << "Running nuebar interactions only " << std::endl;
	nuTypeJnuBeam = 4;
      } else {
	std::cerr << "Invalid Neutrino Type " << nuTypePDG << std::endl;
	std::cerr << "Setting nuType = 0 and generating all species " << std::endl;
	nuTypeJnuBeam = 0;
      }
      ++argv; --argc;
      break;
    case 'r':
      randomSeed = atoi(argv[2]);
      std::cout << std::endl << "Using TRandom3 seed = " << randomSeed << std::endl;
      ++argv; --argc;
      break;
    case 's':
      jnuStem = argv[2];
      ++argv; --argc;
      if(argv[2][0] != '-'){
        jnuStart    = atoi(argv[2]);
        ++argv; --argc;
      } else {
        std::cerr << "Need to specify starting index of jnubeam files " << std::endl;  
        return -1;
      }
      if(argv[2][0] != '-'){
        jnuEnd    = atoi(argv[2]);
        ++argv; --argc;
      } else {
        std::cerr << "Need to specify ending index of jnubeam files " << std::endl;  
        return -1;
      }
      singleFile = false;
      break;
    case 'q':
      rndJnuOrder = atoi(argv[2]);
      ++argv; --argc;
      break;
    case 'Q':
      rndEntryStart = atoi(argv[2]);
      ++argv; --argc;
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

