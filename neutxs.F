************************************************************************
*     -----------------
      subroutine neutxs(ipar,energy,iz,ia,mode,totxs)
*     -----------------
*
*     (Purpose)
*        Return total cross section for a given
*        neutrino type on a specified nuclear target
*
*     (Input)
*        IPAR   : Neutrino type
*                   12  nu_e
*                  -12  nu_e_bar
*                   14  nu_mu
*                  -14  nu_mu_bar
*                   16  nu_tau
*                  -16  nu_tau
*        ENERGY : Neutrino Energy in GeV
*        IZ     : Atomic number
*        IA     : Mass number
*        MODE   : Neut mode
*
*     (Output)
*       TOTXS : Total cross section ( 10^-38cm^2 )
*
*     (Creation Date and Author)
*       2008.11.03 ; R.Tacik
*       2010.11.03 ; P. de Perio - modified to select mode
*
************************************************************************

      IMPLICIT NONE

#include "necardev.h"
#include "necard.h"

      integer ipar, iz, ia, norm, mode, i
      real    energy, totxs, fntotpau
      logical first

      common /CHKINIT/ first
      data first/.true./

C--Read card files once
      if (first) then
         CALL NECARD
         CALL NECARDEV
         first = .false.
      end if
C--

C--   Even though the neut.card file is READ above,
C--   the mode from the call to neutxs subroutine is used
C--   and the mode (and neutrino type) from the
C--   neut.card are NOT used!!
C--   The "average" xsec/nucleon is returned below and then
C--   scaled by number of nucleons in the target.
      CALL NECRSMODE(ipar,mode)


      If (Iz.Eq.1) Then
         Numbndp = 0
         Numbndn = 0
         Numfrep = 1
         Numatom = 1
      Else If (Iz.Gt.1) Then
         Numbndp = Iz
         Numbndn = Ia - Iz
         Numfrep = 0
         Numatom = Ia
         Call Nesetfgparams
      Else
         Write (6,*) "%-Neutxs Called With Negative Value Of Iz=",Iz
         Stop
      End If

      norm = numbndn + numbndp + numfrep

      totxs = float(norm) * fntotpau(ipar,energy)

      return

      END

