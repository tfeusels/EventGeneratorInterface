#include "TObjArray.h"
#include "TGeoMaterial.h"
#include "TGeoNode.h"

#include "TNuTrajectory.h"

// Default Constructor
TNuTrajectory::TNuTrajectory()
{
  fMatList     = 0;
  fNodeList    = 0;
  fNucList     = 0; 

  fMatList     = new TObjArray(0,0);
  fNodeList    = new TObjArray(0,0);
  fNucList     = new TObjArray(0,0);

  fPathMatList = new TObjArray(0,0);

  fXsecFunc    = 0;

  fPathLengths.clear();
  fXsecs.clear();
  fScatCoeffs.clear();

  fNuType      = -1;
  fNuEnergy    = -1;
  fSumIntProb  = -1;
}

// Destructor
TNuTrajectory::~TNuTrajectory()
{
  // delete lists
  if (fMatList) {
    fMatList->Delete();
    delete fMatList;
  }

  if (fNucList) {
    fNucList->Delete();
    delete fNucList;
  }

  if (fNodeList) {
    fNodeList->Delete();
    delete fNodeList;
  }

  fPathLengths.clear();
  fXsecs.clear();
  fScatCoeffs.clear();

  fNuType      = -1;
  fNuEnergy    = -1;
  fSumIntProb  = -1;
}

void TNuTrajectory::Clear(){

  //information from global geometry
  if(fMatList) fMatList->Clear();
  if(fNucList) fNucList->Clear();


  //information specific to last path
  if(fNodeList)fNodeList->Clear();
  if(fPathMatList) fPathMatList->Clear();
  fPathLengths.clear();
  fXsecs.clear();
  fScatCoeffs.clear();

  fNuType      = -1;
  fNuEnergy    = -1;
  fSumIntProb  = -1;
  fMaxIntProb  = -1;

}

void TNuTrajectory::ClearPathInfo(){

  //information specific to last path
  if(fNodeList)fNodeList->Clear();
  if(fPathMatList) fPathMatList->Clear();
  fPathLengths.clear();
  fXsecs.clear();
  fScatCoeffs.clear();
  fPathIntProb.clear();

  fNuType      = -1;
  fNuEnergy    = -1;
  fSumIntProb  = -1;

}

void TNuTrajectory::SetRandomSeed(Int_t seed){
  fRandom.SetSeed(seed);
}

// set geo manager
void TNuTrajectory::SetGeoManager(TGeoManager* theGeoManager){
  fGeoManager = theGeoManager;
}

void TNuTrajectory::SetGlobalVol(TGeoVolume* theGlobalVolume){
  fGlobalVol = theGlobalVolume;
}

void TNuTrajectory::SetStartZ(Double_t z){
  fStartZ = z;
}

void TNuTrajectory::SetMaxIntProb(Double_t d){
  fMaxIntProb = d;
}

void TNuTrajectory::SetBeamWeight(Double_t d){
  fBeamWeight = d;
}

// Add a material to the list of materials for this path
void TNuTrajectory::AddPathMaterial(TGeoMaterial* theMaterial)
{
  if (!fPathMatList) fPathMatList = new TObjArray(128);
    fPathMatList->AddLast(theMaterial);

}


// Add a material to the list of materials ever encountered
void TNuTrajectory::AddMaterial(TGeoMaterial* theMaterial)
{
  if (!fMatList) fMatList = new TObjArray(128);
  
  // add only if material is not already there:
  // if list is empty add it:
  if(fMatList->GetEntries() == 0){
    fMatList->AddLast(theMaterial);
  } else {
    fMatList->Sort();
    if(fMatList->BinarySearch(theMaterial) < 0){
      fMatList->AddLast(theMaterial);
    }
  }
 
  fMatList->Sort();
  
}


// Record intersected node
void TNuTrajectory::AddNode(TGeoNode* theNode)
{
  if (!fNodeList) fNodeList = new TObjArray(128);
  fNodeList->AddLast(theNode);
}
 
// Record pathlength through node
void TNuTrajectory::AddPathLength(Double_t thePathlength){
  fPathLengths.push_back(thePathlength);
}

// set the neutrino cross section function
/*void TNuTrajectory::SetXsecFunc(TNeutXsec* theXsecFunc)
void TNuTrajectory::SetXsecFunc(TXsec* theXsecFunc)
{
  fXsecFunc = theXsecFunc;
  }*/

void TNuTrajectory::SetXsecFunc(std::string nuGenerator)
{
  if(nuGenerator == "neut")
    fXsecFunc = new TNeutXsec();
  else if(nuGenerator == "gibuu")
    fXsecFunc = new TGiBUUXsec();
}


// get a list of all nuclei types in the material
// encountered in the path of the neutrino
void TNuTrajectory::ExtractNuclei()
{

  TObjArrayIter theMaterialIter(fMatList);
  TGeoMaterial* theMat;

  while( (theMat = (TGeoMaterial*)theMaterialIter.Next()) ){

    if(theMat->IsMixture()){

      TGeoMixture* theMix = (TGeoMixture*)theMat;

      for(int ie = 0; ie < theMix->GetNelements(); ie++){

        TGeoElement* theNucleus = theMix->GetElement(ie);

        if(fNucList->GetEntries() == 0){
          fNucList->AddLast(theNucleus);
        } else {
          fNucList->Sort();
          if(fNucList->BinarySearch(theNucleus) < 0){
            fNucList->AddLast(theNucleus);
          }
        }

      } // end loop over elements in mixture
    } else {
      
      TGeoElement* theNucleus = theMat->GetElement();

      if(fNucList->GetEntries() == 0){
        fNucList->AddLast(theNucleus);
      } else {
        fNucList->Sort();

        if(fNucList->BinarySearch(theNucleus) < 0){
          fNucList->AddLast(theNucleus);
        }
      }

    }
    
  } // end loop over mixture array

  fNucList->Sort();
}


void TNuTrajectory::SetNuType(Int_t i)
{
  fNuType = i;
}

void TNuTrajectory::SetNuEnergy(Double_t e)
{
  fNuEnergy = e;
}

void TNuTrajectory::CalcXsecs()
{
  TObjArrayIter theNucleiIter(fNucList);
  TGeoElement* theNucleus;

  while( (theNucleus = (TGeoElement*)theNucleiIter.Next()) ){

    Int_t nProton  = theNucleus->Z();
    Int_t nNeutron = TMath::Nint(theNucleus->A()) - nProton;

    Double_t xsec = -1;

    std::cout << fNuType << " " << fNuEnergy << " " << nProton << " " << nNeutron << std::endl;
    if( nProton == 0 && nNeutron == 0){
      xsec = 0;
    } else {

      xsec = fXsecFunc->getTotalXsec(fNuType,           // neutrino species
				     fNuEnergy,         // neutrino energy
				     nProton,           // protons
				     nProton + nNeutron // A
				     );

      std::cout << xsec << std::endl;
    }

    fXsecs.push_back(xsec);
  }
}


void TNuTrajectory::CalcScatCoefficients(){

  // calculate scattering coefficient (reciprocal of Mean Free Path)
  // for each nucleus type according to the number density and xsec
  TObjArrayIter theMaterialIter(fMatList);
  TGeoMaterial* theMat;

  while( (theMat = (TGeoMaterial*)theMaterialIter.Next()) ){

    Double_t alpha = 0;
    Double_t rhogcm3  = theMat->GetDensity()/(CLHEP::gram/CLHEP::cm3); 
    
    if(rhogcm3 > 0){
      if(theMat->IsMixture()){
	TGeoMixture* theMix = (TGeoMixture*)theMat;
	Double_t*    wfracs   = theMix->GetWmixt();
	
	for(int ie = 0; ie < theMix->GetNelements(); ie++){
	  
	  TGeoElement* theNucleus = theMix->GetElement(ie);
	  Double_t xsec     = GetXsec(theNucleus);
	  Double_t wfrac    = wfracs[ie];
	  Double_t nDens    = rhogcm3 * wfrac * TMath::Na()/theNucleus->A();
	
	  // cross sections are reported in units of 10^{getExponentCm2)
	  // usually, 10^{-38}
	  alpha += xsec *  pow(10,fXsecFunc->getExponentCm2()) * nDens;

	} // end loop over elements in mixture
	
      } else {
	
	
	TGeoElement* theNucleus = theMat->GetElement();
	Double_t xsec     = GetXsec(theNucleus);
	Double_t nDens    = rhogcm3 * TMath::Na()/theNucleus->A();
	
	alpha += xsec * nDens;

      }
    } else {
      alpha = 0;
    }

    fScatCoeffs.push_back(alpha);
     
  } // end loop over material/mixture array


}




// retrieve encountered material information (ever encountered)
TGeoMaterial* TNuTrajectory::GetMaterial(Int_t i)
{
  return (TGeoMaterial*)fMatList->At(i);
}

// retrieve encountered material information (for this path)
TGeoMaterial* TNuTrajectory::GetPathMaterial(Int_t i)
{
  return (TGeoMaterial*)fPathMatList->At(i);
}

TGeoElement* TNuTrajectory::GetNucleus(Int_t i)
{
  return (TGeoElement*)fNucList->At(i);
}

// retrieve encountered node information
TGeoNode* TNuTrajectory::GetNode(Int_t i)
{
  return (TGeoNode*)fNodeList->At(i);
}

// retrieve pathelength within each node
Double_t TNuTrajectory::GetPathLength(Int_t i)
{
  return fPathLengths[i];
}

// retrieve interaction probability on each path
Double_t TNuTrajectory::GetPathIntProb(Int_t i)
{
  return fPathIntProb[i];
}

// retrieve maximum interaction probability based on random sampling
Double_t TNuTrajectory::GetMaxIntProb()
{
  return fMaxIntProb * pow(10,fXsecFunc->getExponentCm2()) ;
}

// retrieve scat coefficient for a material
Double_t TNuTrajectory::GetScatCoefficient(TGeoMaterial* theMat)
{
  Int_t    index = fMatList->BinarySearch(theMat);
  Double_t scatcoeff = -1;

  if(index < 0){
    std::cout << "Error: Material not found" << std::endl;
  } else {
    scatcoeff = fScatCoeffs[index];

  }

  return scatcoeff;

}

// retrieve cross section for a nucleus
Double_t TNuTrajectory::GetXsec(TGeoElement* theNucleus)
{

  Int_t index = fNucList->BinarySearch(theNucleus);
  Double_t xsec = -1;


  
  if(index < 0){
    std::cerr << "Error: Nucleus not found" << std::endl;
  } else {
    xsec = fXsecs[index];
  }

  
  return xsec;
  
}

// retrieve number of materials;
Int_t TNuTrajectory::GetNMaterials()
{
  return fMatList->GetEntries();
}

// retrieve number of nuclei;
Int_t TNuTrajectory::GetNNuclei()
{
  return fNucList->GetEntries();
}

// retrieve number of nodes;
Int_t TNuTrajectory::GetNNodes()
{
  return fNodeList->GetEntries();
}  


// retrieve number of nodes;
void TNuTrajectory::Swim(TNuVector* theNuVector)
{
  
  ClearPathInfo();

  fDir   = theNuVector->getNuDir().Unit();

  Double_t displacement = (theNuVector->getNuPos().Z() - fStartZ)/fDir.CosTheta();

  this->SetNuType(theNuVector->getNuType());
  this->SetNuEnergy(theNuVector->getNuEnergy());
 
  fStartPoint = theNuVector->getNuPos() - fDir*displacement;

  //  Double_t globalStart[3] = {}
  //  Double_t localStart[3]  = {0,0,0};
  //  TranslateGlobalToLocal(fTopNode, fTopVolume);

  Double_t spoint[3] = {fStartPoint.X(), fStartPoint.Y(), fStartPoint.Z()};
  Double_t sdir[3]   = {fDir.X(), fDir.Y(), fDir.Z()};

  fGeoManager->InitTrack(spoint[0], spoint[1], spoint[2],
			 sdir[0]  , sdir[1],   sdir[2]   );

  TVector3 lastPoint; lastPoint.SetXYZ(spoint[0], spoint[1], spoint[2]);
  Double_t point[3] = {0,0,0};
  Bool_t isGeom = true, cross = true;

  // set the beam weight
  SetBeamWeight((Double_t)theNuVector->getNuWeight());

  do {
    // get the current node and volume
    TGeoNode     *cnode = fGeoManager->GetCurrentNode(); 
    TGeoVolume   *cvol  = fGeoManager->GetCurrentVolume(); 
    TGeoMaterial *cmat  = cvol->GetMedium()->GetMaterial();

    // add current node/material information to list
    this->AddNode(cnode);
    this->AddMaterial(cmat);
    this->AddPathMaterial(cmat);

   // calculate step size
    TGeoNode* theNextNode = fGeoManager->FindNextBoundary(); 

    // make the step
    theNextNode = fGeoManager->Step(isGeom, cross);

    const Double_t* tmppoint = fGeoManager->GetCurrentPoint();
    point[0] = tmppoint[0];
    point[1] = tmppoint[1];
    point[2] = tmppoint[2];

    TVector3 currentPoint; 
    currentPoint.SetXYZ(point[0], point[1], point[2]);

    this->AddPathLength((lastPoint-currentPoint).Mag());

    lastPoint = currentPoint;



  } while( fGlobalVol->Contains(point) );

  return;
}  

Double_t TNuTrajectory::GetSumIntProb(){
  return fSumIntProb;
}

void TNuTrajectory::CalcSumIntProb(){

  
  fSumIntProb = 0;
  fPathIntProb.clear();

  for(int ip = 0; ip < (Int_t)fPathLengths.size(); ip++){
    Double_t      pathLength  = this->GetPathLength(ip);
    TGeoMaterial* theMaterial = this->GetPathMaterial(ip);

    Double_t scatCoeff        = 0;

    scatCoeff = this->GetScatCoefficient(theMaterial);

    // scat coeff is in cm-1 so path length must be in cm
    Double_t pathIntProb      = scatCoeff * pathLength/CLHEP::cm;
    
    fPathIntProb.push_back(pathIntProb);
    fSumIntProb  += pathIntProb;
  }
  
  return;

}

// Calculate raw interaction probability per 1e21 POT
Double_t TNuTrajectory::GetRawProbability(){
  if(fSumIntProb == -1){
    std::cout << "Error: have not yet calculated sum interaction probability"
	      << std::endl;
  }

  Double_t intProb = fSumIntProb * pow(10,fXsecFunc->getExponentCm2()) *fBeamWeight;
  return intProb;
}

// determine if interaction actually happens:
Bool_t TNuTrajectory::DrawInteraction(){
 
  Double_t relativeIntProb = GetRawProbability()/GetMaxIntProb();
 
  if(relativeIntProb > 1)
    std::cout << "Error: interaction probability is greater than max" 
	      << std::endl;
  return (fRandom.Rndm() < relativeIntProb);
  
}

TLorentzVector TNuTrajectory::DrawVertexFromPath(TGeoElement*& theTarget){

  if(fSumIntProb == -1){
    std::cout << "Error: have not yet calculated sum interaction probability"
	      << std::endl;
  }

  if(fMaxIntProb < 0){
    std::cout << "Error: have not yet calculated max interaction probability"
	      << std::endl;
  }

  Double_t randomDraw = fSumIntProb * fRandom.Rndm();

  Double_t intProb     = 0;
  Int_t    iPath       = 0;
  TVector3 theVertex   = fStartPoint;
  Double_t theTime     = 0;
  TGeoMaterial* theMat = 0;
  Double_t c           = TMath::C()*CLHEP::m/CLHEP::s; // TMath::C() gives c in m/ssec
  
  while( (intProb + this->GetPathIntProb(iPath)) < randomDraw){
    intProb   += this->GetPathIntProb(iPath);
    theVertex += fDir * this->GetPathLength(iPath);
    theTime   += c*this->GetPathLength(iPath);
    iPath++;
  } 

  if(iPath > fPathIntProb.size()){
    std::cout << "Error: selected path is out of range " << iPath << std::endl;
  }
  
  Double_t remainingProb = (randomDraw-intProb)/this->GetPathIntProb(iPath);
  theVertex += fDir * remainingProb * this->GetPathLength(iPath);
  theTime   += c    * remainingProb*this->GetPathLength(iPath);

  theMat     = this->GetPathMaterial(iPath);

  if(theMat->IsMixture()){
    TGeoMixture* theMix = (TGeoMixture*)theMat;
    Double_t*    wfracs   = theMix->GetWmixt();
    Double_t rhogcm3  = theMat->GetDensity()/(CLHEP::gram/CLHEP::cm3); 

    Double_t theDraw = fRandom.Rndm() * this->GetScatCoefficient(theMat);
    Double_t sumAlpha=0;
    Int_t    ie = 0;


    do {
	  
      TGeoElement* theNucleus = theMix->GetElement(ie);
      Double_t xsec     = GetXsec(theNucleus) * pow(10,fXsecFunc->getExponentCm2());
      Double_t wfrac    = wfracs[ie];
      Double_t nDens    = rhogcm3 * wfrac * TMath::Na()/theNucleus->A();


      sumAlpha += nDens * xsec;
      theTarget = theNucleus;

      ie++;

    } while (sumAlpha < theDraw);// end loop over elements in mixture

  } else {
    theTarget = theMat->GetElement();
  }
  TLorentzVector theFourVertex(theVertex, theTime);
  return theFourVertex;

}

TNuVertex* TNuTrajectory::DrawVertexFromPath(){

  TGeoElement*   theTargetNucleus = 0;
  TLorentzVector theVertex   = DrawVertexFromPath(theTargetNucleus);
  TNuVertex*     theNuVertex = new TNuVertex(theVertex, theTargetNucleus);

  return theNuVertex;
  

}
