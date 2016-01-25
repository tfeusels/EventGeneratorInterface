#ifndef XSEC_H
#define XSEC_H

/**************************************************************************************
 * TXsec.h
 *
 * Abstract class which will work for ANY generator and makes calls in event_rate and
 * getev generic!
 * 
 *
 * Calls the total cross section (implemented in derived classes for each generator)
 * Control for obtaining cross sections for exclusive channels is through
 * the card files (neut.card). If exclusive channels are specified there,
 * getTotalXsec will return the cross section for the exclusive process
 *
 *
 **************************************************************************************/

class TXsec {

 public:
  
  TXsec(){};

  ~TXsec(){ };

  // Common getters and setters for all generators
  // Initial value will be set for each generator during construction!
  void  setExponentCm2(Int_t i) { fExponentCm2 = i;    };
  Int_t getExponentCm2()        { return fExponentCm2; };

  /* The function expects:
   * 1. neutrino species +/- 12/14/16 for nue/numu/nutau
   * 2. energy in CLHEP units (MeV)
   * 3. Z/A of target nucleus
   */
  virtual float getTotalXsec(int nutype, float nu_energy, int z_target, int a_target) = 0; 

 protected:

  // The variable fExponentCm2 stores the log_10 of the units of the
  // returned cross section in cm2, i.e. -38 means that the cros sections are
  // returned in units of 10^{-38} cm^2.
  int fExponentCm2;

};

#endif
