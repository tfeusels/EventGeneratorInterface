//////////////////////////////////////////////////////////////////////////
// 
// TNeutOutput.h
// Root output class; takes generated information and shoves it into:
// Neut output format (fFormat = 0)
// RooTracker format  (fFormat = 1)
//
//////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TBits.h"
#include "TObjString.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TNuVertex.h"
#include "TNuVector.h"
#include "TNeutEventGenerator.h"
#include <iostream>

class TNeutOutput {

public:
  TNeutOutput(){};
  TNeutOutput(TString filename, Int_t format=0, Float_t fluxVersion=0);
  
  void OpenFile(TString filename);
  void InitTree(Float_t fluxVersion);
  void SetOutputFormat(Int_t i) { fFormat = i; } // 0 for Neut format
                                                 // 1 for RooTracker format
  void FillTree(Int_t                eventNumber,
		TNuVector*           theNuVector,
		TNuVertex*           theNuVertex,
		//		TNeutEventGenerator* theNuFinalState);
		TNuEventGenerator* theNuFinalState);

  void ClearNeutTree();
  void ClearRooTrackerTree();
  void FillNeutTree(Int_t                      eventNumber,
		    TNuVector*                 theNuVector,
		    TNuVertex*                 theNuVertex,
		    TNeutEventGenerator*       theNuFinalState);
  
  void FillRooTrackerTree(Int_t                eventNumber,
			  TNuVector*           theNuVector,
			  TNuVertex*           theNuVertex,
			  //			  TNeutEventGenerator* theNuFinalState);
			  TNuEventGenerator* theNuFinalState);

  void SetTreeWeight(double potToGenerate=1);
  
  void WriteTree();
  void CloseFile();
 
private:
  TFile*       fOutputFile;
  TTree*       fOutputTree;
  Int_t        fFormat;
  TDatabasePDG fPDGDatabase;
  
  // Neut Format
  /// The event number.
  int   fNev;
  
  /// The event vertex
  float fPos[3];
  
  /// The interaction mode of the event.
  int   fMode;
  
  /// The number of particles in the neut kinematics list.
  int   fNumnu;
  
  /// The PDG code for the particles in the kinematics list.
  int   fIpnu[100];   //[Numnu]
  
  /// The absolute value of the momentum (GeV/c)
  float fAbspnu[100];   //[Numnu]
  
  /// The particle momentum. (GeV/c)
  float fPnu[100][3];   //[Numnu]
  
  /// The number of particles leaving the nucleus.  This is the total number
  /// of particles that should be tracked by G4.
  int   fNpar;
  
  /// The particle ID of particles leaving the nucleus.  These particles
  /// should be tracked by G4.
  int   fIpv[100];   //[Npar]
  
  /// Flag if this particle should be tracked in the detector.  Track if
  /// this is not zero.
  int   fIcrnv[100];   //[Npar]
  
  /// The three momentum of particles leaving the nucleus.
  float fPmomv[100][3];   //[Npar]
  
  /// The G3 particle code for the neutrino parent particle.
  int   fppid;

  /// The decay position of the neutrino parent from the beam line
  /// simulation.  This is the position with respect to the target.
  float fxpi[3];
  
  /// The direction of the neutrino parent from the beam line simulation.
  float fnpi[3];
  
  /// The cosine of the angle between the neutrino parent direction and the
  /// beam direction
  float fcospibm;
  
  /// The momentum of the neutrino parent.  This is in units of GeV/c
  float fppi;
  

  
  /// The momentum of the neutrino parent particle at the production point.
  float fppi0;

  /// The production position for neutrino's parent particle from the beam
  /// line simulation.
  float fxpi0[3];
  
  /// The production direction of the neutrino parent from the beam line
  /// simulation.
  float fnpi0[3];
  
  /// The cosine of the angle between the neutrino parent production
  /// direction and the beam direction.  This is the direction of the
  /// particle at the target.
  float fcospi0bm;
  

  
  
  //
  // RooTrackerFormat
  //
  
  /// The generator-specific event flags.
  TBits*       fEvtFlags;
  
  /// The generator-specific string with the 'event code'
  TObjString*  fEvtCode;
  
  /// The sequence number of the event (the event number).
  int         fEvtNum;

  /// Native NEUT flag for bound or unbound proton
  int         fNEibound;
  
  /// The cross section for the event (1E-38 cm2)
  double      fEvtXSec;
  
  /// The differential cross section for the event kinematics 
  /// (1E-38 cm2/{K^n})
  double      fEvtDXSec;
  
  /// The weight for the event
  double      fEvtWght;
  
  /// The probability for the event (given the cross section, path lengths,
  /// etc.).
  double      fEvtProb;
  
  /// The event vertex position in detector coord syst (in meters and
  /// seconds).
  double      fEvtVtx[4];
  
  /// The number of particles in the particle arrays to track
  int         fStdHepN;
  
  /// The maximum number of particles that can be in the particle arrays.
  static const int kNPmax = 4000;
  
  /// The PDG codes for the particles to track.  This may include generator
  /// specific codes for pseudo particles.
  int         fStdHepPdg[kNPmax]; //[fStdHepN]
  
  /// The a generator specific status for each particle.  Particles with a
  /// status equal to 1 will be tracked.
  ///
  /// The HEPEVT status values are as follows:
  /// - 0 : null entry.
  /// - 1 : an existing entry, which has not decayed or fragmented. This is
  ///    the main class of entries, which represents the `final state' given
  ///    by the generator.
  /// - 2 : an entry which has decayed or fragmented and is therefore not
  ///    appearing in the final state, but is retained for event history
  ///    information.
  /// - 3 : a documentation line, defined separately from the event
  ///    history. This could include the two incoming reacting particles,
  ///    etc.
  /// - 4 to 10 :
  ///    undefined, but reserved for future standards. 
  /// - 11 to 200 : at the disposal of each model builder for constructs
  ///    specific to his program, but equivalent to a null line in the
  ///    context of any other program.
  /// - 201 and above : at the disposal of users, in particular for event
  /// tracking in the detector.
  ///
  /// The GENIE generator approximately follows the HEPEVT status codes.
  /// As of July 2008, the status values found the GENIE source code are:
  ///   - -1 -- Undefined particle
  ///   -  0 -- An initial state particle.
  ///   -  1 -- A stable final state particle to be tracked.
  ///   -  2 -- An intermediate particle that should not be tracked.
  ///   -  3 -- A particle which decayed and should not be tracked.  If 
  ///            this particle produced daughters to be tracked, they will
  ///            have a state of 1.
  int         fStdHepStatus[kNPmax]; //[fStdHepN]
  
  /// The position (x, y, z, t) (fm, second) of the particle in the nuclear
  /// frame
  double      fStdHepX4[kNPmax][4]; //[fStdHepN]
  
  /// The 4-momentum (px, py, pz, E) of the particle in the LAB frame (GeV)
  double      fStdHepP4[kNPmax][4]; //[fStdHepN]
  
  /// The particle polarization vector.
  double      fStdHepPolz  [kNPmax][3]; //[fStdHepN]
  
  /// The index of the first daughter of the particle in the arrays.
  int         fStdHepFd[kNPmax]; //[fStdHepN]
  
  /// The index last daughter of the particle in the arrays.
  int         fStdHepLd[kNPmax]; //[fStdHepN]
  
  /// The index of the first mother of the particle in there arrays.
  int         fStdHepFm[kNPmax]; //[fStdHepN]
  
  /// The index of the last mother of the particle in the arrays.
  int         fStdHepLm[kNPmax]; //[fStdHepN]
  



  // NEUT native VCWORK information (Added Sep. 24, 2010)
  
  int fNEneutmode;

  /// The maximum number of particles that can be in the particle arrays.
  static const int kNEmaxvc = 100;
  
  int    fNEnvc;
  //float  fNEposvc[3];
  int    fNEipvc[kNEmaxvc];
  //float  fNEamasvc[kNEmaxvc];
  float  fNEpvc[kNEmaxvc][3];
  int    fNEiorgvc[kNEmaxvc];
  int    fNEiflgvc[kNEmaxvc];
  int    fNEicrnvc[kNEmaxvc];
  //float  fNEtimvc[kNEmaxvc];
  //float  fNEposivc[kNEmaxvc][3];
  //int    fNEivtivc[kNEmaxvc];
  //float  fNEposfvc[kNEmaxvc][3];
  //int    fNEivtfvc[kNEmaxvc];

  // Additional cross section kinematic info
  float  fNEcrsx   ;
  float  fNEcrsy   ;
  float  fNEcrsz   ;
  float  fNEcrsphi ;

  // NEUT FSIHIST pion interaction history

  static const int kNEmaxvert = 100;
  int    fNEnvert;
  float  fNEposvert[kNEmaxvert][3];
  int    fNEiflgvert[kNEmaxvert];

  static const int kNEmaxvertp = 300;
  int    fNEnvcvert;
  float  fNEdirvert[kNEmaxvertp][3];
  float  fNEabspvert[kNEmaxvertp];
  float  fNEabstpvert[kNEmaxvertp];
  int    fNEipvert[kNEmaxvertp];
  int    fNEiverti[kNEmaxvertp];
  int    fNEivertf[kNEmaxvertp];


  // FLUX INFORMATION


  /// The PDG code of the particle which created the parent neutrino.
  int         fNuParentPdg;
  
  /// The interaction mode at the vertex which created the parent neutrino.
  /// This is normally the decay mode of the parent particle.
  int         fNuParentDecMode;
  
  /// The 4 momentum of the particle at the vertex which created the parent
  /// neutrino.  This is normally the momentum of the parent particle at the
  /// decay point.
  double      fNuParentDecP4[4];
  
  /// The position (meters, seconds) of the vertex at which the neutrino was
  /// created.  This uses the target as the origin.
  double      fNuParentDecX4[4];
  
  float      fNuCospibm;

  float fNuNorm;

  /// The momentum (GeV) of the parent particle at it's production point.
  double      fNuParentProP4[4];

  /// The position (meters, seconds) of the parent particle at it's
  /// production point.  This uses the target as the origin.
  double      fNuParentProX4[4];
  
  float      fNuCospi0bm;

  /// The vertex ID of the parent particle vertex.
  int         fNuParentProNVtx;

  // The neutrino position at flux plane
  float      fNuRnu, fNuXnu[2]; 

  // Detector ID
  int fNuIdfd;


  // primary particle information
  // primary particle ID
  int         fNuGipart;
  // primary particle starting point
  float      fNuGpos0[3];
  // primary particle direction at starting point
  float      fNuGvec0[3];
  // momentum of the primary particle at the starting point
  float      fNuGamom0;

  
  // Flux Secondary Interation History Information (>=2011a)

  // Maximum number of interaction steps
  static const int kNgmax = 12; 

  // Actual number of interaction steps
  int   fNuNg;
  // Momentum of intermediate particle
  float fNuGp[kNgmax][3];
  float fNuGcosbm[kNgmax];
 // Position of intermediate particle
  float fNuGv[kNgmax][3];
  // PDG particle ID 
  int fNuGpid[kNgmax];
  // GEANT Interaction mechanism code
  int fNuGmec[kNgmax]; 
 
  // Normalization and Transfer Information
  float fNuEnusk, fNuNormsk, fNuAnorm;

  int    fNuGmat[kNgmax]; 
  float  fNuGdistc[kNgmax];  
  float  fNuGdistal[kNgmax]; 
  float  fNuGdistti[kNgmax]; 
  float  fNuGdistfe[kNgmax]; 




  // Beam parameter information (>=2010d)
  
  // jnubeam flux version
  float fNuVersion;

  // beam tune ID #
  int fNuTuneid;

  int fNuNtrig;
  
  // Interaction model ID
  int fNuPint;

  // Beam center position and angle (divergence?)
  float fNuBpos[2], fNuBtilt[2];
  
  // Beam RMS width
  float fNuBrms[2];

  // Beam optics parameters
  float fNuEmit[2], fNuAlpha[2];

  // Horn currents
  float fNuHcur[3];

  int fNuRand;
  int fNuRseed[2];

  // File information
  long fNuFluxEntry;
  TObjString* fNuFileName;
};
