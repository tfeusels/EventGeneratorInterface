!************************************************************************
!*     -----------------
subroutine gibuuxs(nutype,energy,target_z,target_a,mode,totxs)
!*     -----------------
!*
!*     (Purpose)
!*        Return total cross section for a given
!*        neutrino type on a specified nuclear target
!*        and mode
!*
!*     (Input)
!*        NUTYPE   : Neutrino type
!*                   12  nu_e
!*                  -12  nu_e_bar
!*                   14  nu_mu
!*                  -14  nu_mu_bar
!*                   16  nu_tau
!*                  -16  nu_tau
!*        ENERGY : Neutrino Energy in GeV
!*        TARGET_Z     : Atomic number
!*        TARGET_A     : Mass number
!*        MODE   : Gibuu mode (like NEUT, also overwrite any specification in .job)
!*                             --> .job file is for model/generator settings !!!
!*
!*     (Output)
!*       TOTXS : Total cross section ( 10^-38cm^2 )
!*
!*     (Creation Date and Author)
!*       2014.03.25 ; T.Feusels
!*
!************************************************************************
!module gibuuXsCalc 

  !BASED on init_neutrino.f90 and GiBUU.f90
  use particleDefinition
  use nucleusDefinition, only: tnucleus
  use eN_eventDefinition
  use eN_event
  use neutrino_IDTable
  use idtable
  use neutrinoXsection
  use neutrino_IDTable
  use random, only : rn
  use nucleusDefinition
  use NuclearPDF, only : SetNuclearPDFA
  use neutrinoInfoStorage
  use neutrinoProdInfo, only: neutrinoProdInfo_Init,neutrinoProdInfo_Store
  use inputGeneral , only : fullEnsemble
  use pauliBlockingModule, only : checkPauli
  use nucleus, only: getTarget, initNucleus  !! hack needed: nucleus.f90 L20: made InitNucleus public
  use inputGeneral,  only: eventType, readinputGeneral, length_perturbative, length_real, numEnsembles, printParticleVectors
  use baryonPotentialModule, only: HandPotentialToDensityStatic
  use particleProperties, only: initParticleProperties
  use checks, only: ChecksSetDefaulSwitches

  implicit none

  !! I/O
  integer, intent(in) :: nutype
  integer, intent(in) :: target_z
  integer, intent(in) :: target_a
  integer, intent(in) :: mode
  real, intent(out) :: totxs
  real, intent(in) :: energy

  logical :: debugFlag=.false.
  logical :: initFlag = .true.

  !! target and initial state particles
  type(particle), Allocatable :: realParticles(:,:)
  type(particle), Allocatable :: pertParticles(:,:)
  type(tnucleus), pointer :: targetNuc => NULL()
  integer :: lengthReal = 0  ! maximum number of real particles per ensemble
  integer :: lengthPert = 0  ! maximum number of perturbative particles per ensemble

  !! for xsec calculation
  type(electronNucleon_event) :: eNev0, eNev1
  integer, parameter :: max_finalstate_ID=37
  type(electronNucleon_event), dimension(1:max_finalstate_ID) :: eNev
  
  !! output kinematics from xsec
  type(particle),dimension(1:max_finalstate_ID,1:2),target :: OutPart
  type(particle),dimension(1:20),target :: OutPartDIS
  type(particle),dimension(1:3), target :: OutPart2pi
  type(particle),dimension(:),pointer :: pOutPart

  !! total cross section
  real :: sigtot
  real*8, dimension(1:max_finalstate_ID) ::  sigma=0.
  real,dimension(0:max_finalstate_ID)    :: sigabsArr

  integer :: i, j, k
  integer :: flavor_ID
  integer :: process_ID=2
  logical :: MCmode=.false.
  logical :: notPauliBlocked
  integer :: whichReal_nucleon, numNucleons
  integer :: numtry=0
  integer :: nuXsectionMode=0
  real*8 :: raiseVal
  logical, dimension(1:max_finalstate_ID) :: includeK   !for switching included channels in xsec on and off

  character*(*), dimension(3), parameter :: sProcess = (/"EM","CC","NC"/)
  character*(*), dimension(3), parameter :: sFamily  = (/"e  ","mu ","tau"/)
  logical, save :: includeQE=.true.
  logical, save :: includeDELTA=.true.
  logical, save :: includeRES=.true.
  logical, save :: include1pi=.false.
  logical, save :: include2pi=.false.
  logical, save :: includeDIS=.false.
  logical, save :: include2p2hQE=.false.
  logical, save :: include2p2hDelta=.false.
  real, save :: Enu_upper_cut=200.
  real, save :: Enu_lower_cut=0.
  integer, parameter :: max_Hist = 8 ! number of different histograms
  logical, save :: realRun=.false.
  logical, dimension(0:max_Hist), save :: includeHist
  real, save :: sigmacut=10e-4
  logical, save :: readinputflag = .true.
  integer, save :: nuExp=0  


  ! DEAL WITH READING INPUT/CARD file, settings and channel choice
  call readInputGeneral             !read input namelist of .job file
  call initParticleProperties       !define variables baryon and meson
  call ChecksSetDefaulSwitches(eventType)

  !! fullensemble (num_ensemble test particles per nucleon for BUU calculations
  !! See GIBUU manual.
  !! Q: Do I also need the ensemble for the cross section calculation?
  !! A: If the test particles are different (well, different density/position)
  !!   => BUT does that affect the xsec? Well, try out and compare with different
  !!      ensembles. Will be MUCH faster if num_ensemble is one here, but necessary
  !!      for initial state kinematics for later BUU calculations anyway.


  !! set up target according to .job file (eg. fermiMotion, density, fermiMomentum, Binding)
  targetNuc => getTarget()          !set up target resting at 0 (and read target namelist from .job)
  !! Overwrite A and Z with input from subroutine!
  targetNuc%mass = target_a
  targetNuc%charge = target_z
  call initNucleus(targetNuc,targetNuc%fermiMomentum_input)  !! recalc stuff based on new A and Z
  !! TODO : shouldn't I just rewrite getTarget?
  


  lengthReal = targetNuc%mass ! Real particles fixed by nucleons in target
  
  if (targetNuc%ReAdjustForConstBinding) then
     write(*,*) 'we now have to readjust the density'
     call HandPotentialToDensityStatic(targetNuc)
  end if
  
  !! realRun = false ( default behaviour: initialize initial state particles from nu-nucl interaction
  !!                                      as perturbative particles for the BUU equations to propagate)
  lengthPert =  max(100,10*targetNuc%mass)
  
  !!

  if (length_perturbative >= 0) lengthPert = length_perturbative
  if (length_real >= 0)         lengthReal = length_real
  
  !...Allocate the vectors
  Allocate(realparticles(1:numEnsembles,1:lengthReal))
  Allocate(pertparticles(1:numEnsembles,1:lengthPert))

  !! fill the nensemble x nNucleon matrix
  call SetUpTarget()
  
  !!!!!!!!!!!!!
  !call init_Neutrino(realParticles: NUCLEUS with NUCLEONS with NENSEMBLE/NUCLEON,
  !                   pertParticles: INITIAL STATE particles given by xsec for channel k,
  !                   energyRaiseFlag: for energy scan => FALSE here,num_runs_sameEnergy = 1 HERE,
  !                   targetNuc)
  !!!!!!!!!!!!!
  if(initFlag) then
     call readInput
     call DoInit
     initFlag=.false.
  end if

  sigabsArr = 0.

  !loop to determine numtry (number of test particles)
  ! TF: different from GiBUU initNeutrino: I need totxs/nucleus => average over ensemble, not over nucleons
  !numtry=numEnsembles

  
  !! call neutrinoProdInfo_Init(numtry) !! initializes storage of neutrino event (not needed for xsec
  call SetNuclearPDFA(targetNuc%mass) 

  ! set the overall kinematics (most of it as dummy):
  flavor_ID = (abs(nutype)-10)/2       !! PDG to GiBUU convention
  call eNev_SetProcess(eNev0, sign(process_ID,nutype),flavor_ID)   !! QUESTION/TODO: loop over NC/CC? YES!!!
  
  loopOverEnsemble: do i = lBound(realParticles,dim=1),uBound(realParticles,dim=1)
     !print out status info
     if(mod(i,50)==0.) write(*,*) 'now starting ensemble ',i
     
     loopVector: do j = lBound(realParticles,dim=2),uBound(realParticles,dim=2)

        if(realParticles(i,j)%ID.ne.nucleon) cycle
        !if (energy.lt.Enu_lower_cut .or. energy.gt.Enu_upper_cut) then
        !        countcutoff=countcutoff+1
        !        cycle
        !end if
        
        eNev1 = eNev0
        call eNev_init_nuStep1(eNev1,realParticles(i,j))
        write(*,*) 'eNev ', eNev1%nucleon
        !write(*,*) 'eNev2 ', eNev1%nucleon2


        sigma=0.
        if (MCmode) call SetXsecMC(eNev1,dfloat(energy),nuXsectionMode)
        
        call resetNumberGuess()
     
        !calculate cross sections
        particle_ID_loop: do k=1, max_finalstate_ID
           
           if (.not.includeK(k)) cycle
           
           eNev(k) = eNev1
           
           select case (k)
           case(DIS_CH)
              pOutPart => OutPartDIS(:)
           case(twoPion)
              pOutPart => OutPart2pi(:)
           case DEFAULT
              pOutPart => OutPart(k,:)
           end select
           
           !! calc total cross section for neutrino and target (stored in eNev) and channel k
           !! also calls readinput to read in the namelists nl_neutrinoxsectinon (and nl_SigmaMC)
           call Xsec_integratedSigma(eNev(k), k, &
                &  .false.,raiseVal,pOutPart,sigma(k), dfloat(energy))  !!convert float to real(8)!

           ! Pauli blocking
           if(sigma(k).gt.0.) then
              call checkPauli(pOutPart,realParticles,notPauliBlocked)
              if(.not.notPauliBlocked) sigma(k)=0.
           end if
           
           write(*,*) 'sig(k) ', k, ' : ', sigma(k)
        end do particle_ID_loop
        
        ! total cross section for all channels for one target nucleon and one ensemble
        sigtot=sum(sigma)
        !if (debugflag) write(*,'(A,g13.5)') ' sigtot= ', sigtot
        write(*,'(A,g13.5)') ' sigtot= ', sigtot
                
        ! average over all nucleons and ensemble
        sigabsArr(0)  = sigabsArr(0)  + sigtot/float(numtry)
        sigabsArr(1:) = sigabsArr(1:) + sigma(1:)/float(numtry)
                
     end do loopVector
  end do loopOverEnsemble

  totxs = sigabsArr(0)   !should be same as sigtot for nensemble == 1 (and also for nensemble == n?)
  write(*,'(A,g13.5)') ' totxs= ', totxs
  
  !! can't include from main program GiBUU.f90
  !! same for readInput and DoInit from InitNeutrino.f90 because they're within the subroutine
  contains
    
    !! copy-paste from GiBUU.f90 (with UpdateEnergies = true as is default for neutrinos)
    !! fills the realParticles with nucleons from targetNuc !!!
    subroutine SetUpTarget()
      use initNucleus_in_PS,only : initNucPhaseSpace
      use densityModule, only: updateDensity
      use coulomb, only: updateCoulomb
      use yukawa, only: updateYukawa
      use propagation, only: updateVelocity
      use energyCalc, only: updateEnergies
      
      call initNucPhaseSpace(realparticles,targetNuc)
      call updateDensity(realParticles)
      call updateCoulomb
      call updateYukawa(.true.)
      call updateVelocity(realParticles)
      call updateEnergies(realParticles)
      
    end subroutine SetUpTarget

    subroutine DoInit()
      !---------------------------------------------------------------------
      ! setting up some arrays for switching on/off the channels
      !---------------------------------------------------------------------
      includeHist = (/.true.,includeQE, includeDELTA, includeRES, &
           include1pi, includeDIS, include2p2hQE, include2p2hDelta, &
           include2pi /)
      
      includeK = .false.
      do k=1,max_finalstate_ID
         
         select case (k) ! === check for inclusion of process
         case (nucleon)
            if (.not.includeQE) cycle
         case (delta)
            if(.not.includeDELTA) cycle
         case (P11_1440:F37_1950)
            if(.not.includeRES) cycle
         case (onePionCH_n,onePionCH_p)
            if(.not.include1pi) cycle
         case (DIS_CH)
            if(.not.includeDIS) cycle
         case (QE2p2h)
            if(.not.include2p2hQE) cycle
         case (delta2p2h)
            if(.not.include2p2hDelta) cycle
         case (twoPion)
            if(.not.include2pi) cycle
         end select
         
         select case (k) ! === check for inclusion of resonance
         case (S11_2090,D13_2080,G17_2190,P11_2100,P13_1900,F15_2000,S31_1900,D33_1940,D35_1930,D35_2350,P31_1750,F35_1750)
            cycle
         end select
         
         includeK(k) = .true.
      end do
      
    end subroutine DoInit
    
    
    !!read from neutrino_induced namelist
    !! I prefer to "use" this subroutine, but has to be moved out of initNeutrino then
    !! this is the slimmed down version without printing
    subroutine readInput
      use output
      
      integer :: ios
      
      NAMELIST /neutrino_induced/  process_ID,flavor_ID,nuXsectionMode,nuExp, &
           & debugFlag,includeQE,includeDELTA,includeRES,include1pi,includeDIS,&
           & include2p2hQE, include2p2hDelta, include2pi, &
           & sigmacut, realRun, &
           & Enu_lower_cut, Enu_upper_cut
      
      
      if (.not.readinputflag) return
      
      call Write_ReadingInput('neutrino_induced',0)
      rewind(5)
      read(5,nml=neutrino_induced,IOSTAT=ios)
      call Write_ReadingInput("neutrino_induced",0,ios)
      
      
      readinputflag = .false.
      
    end subroutine readInput
  
    
    
    

end subroutine gibuuxs
