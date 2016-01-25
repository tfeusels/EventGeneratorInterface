************************************************************************
*     -----------------
      subroutine neutev(ipar,numom,iz,ia,ierr)
*     -----------------
*
*     (Purpose)
*        Generate neutrino interaction for a given
*        neutrino type with a given momentum and
*        direction on a specified nuclear target
*
*     (Input)
*        IPAR     : Neutrino type
*                     12  nu_e
*                    -12  nu_e_bar
*                     14  nu_mu
*                    -14  nu_mu_bar
*                     16  nu_tau
*                    -16  nu_tau
*        NUMOM[3] : Neutrino Energy in GeV
*        IZ       : Atomic number
*        IA       : Mass number
*
*     (Output)
*       via NEWORK common block
*       ierr      : error code (0: ok, 1: error)
*
*     (Creation Date and Author)
*       2008.11.05 ; R.Tacik
*
************************************************************************

      IMPLICIT NONE

#include "necardev.h"
#include "necard.h"
#include "nework.h"
#include "vcwork.h"
#include "fsihist.h"

      integer ipar, iz, ia, ierr
      real    numom(3), pos(3), pneu(3)
      logical first

      common /CHKINIT/ first

C--Read card files once

      if (first) then
         CALL NECARD
         CALL NECARDEV
         first = .false.
      end if

      ierr = 0

C--Set nuclear target parameters
      if (iz.eq.1) then
         NUMBNDP = 0
         NUMBNDN = 0
         NUMFREP = 1
         NUMATOM = 1
      else if (iz.gt.1) then
         NUMBNDP = iz
         NUMBNDN = ia - iz
         NUMFREP = 0
         NUMATOM = ia
         call nesetfgparams
      else
         write (6,*) "%-neutev called with negative value of iz=",iz
         ierr = 1
C         return
         stop
      end if

C--Convert units:
C  numom is in GeV, for compatibility with neutxs
C  pneu, input to nevecgen, wants to be in MeV

      PNEU(1) = numom(1)*1000.
      PNEU(2) = numom(2)*1000.
      PNEU(3) = numom(3)*1000.

C--nevecgen wants a position for the initial vertex
C  set this to zero for now

      POS(1) = 0.
      POS(2) = 0.
      POS(3) = 0.

C--Generate event

      call NEVECGEN(IPAR,POS,PNEU,IERR)

      return

      END
