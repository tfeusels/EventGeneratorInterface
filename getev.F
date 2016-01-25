      implicit none

#include "nework.h"

      integer nu_id/14/
      integer i,j,iz,ia
      real    e(3)

      iz = 8
      ia = 16

      e(1) = 0.
      e(2) = 0.
      e(3) = 0.5


      do i = 1, 10

         call neutev(nu_id,e,iz,ia)

C         write (6,*) 'Event number ',i
C         write (6,*) 'Interaction mode ',modene
C         write (6,*) 'Number of particles ',numne

c         do j = 1, numne
c            write (6,20) j,ipne(j),iorgne(j),icrnne(j),iflgne(j)
c 20         format('Flags for particle ',i2,' : ',2x,i6,3(2x,i3))
c            write (6,30) pne(1,j),pne(2,j),pne(3,j)
c 30         format('Momentum x, y, z : ',3(2x,e12.5))
c         end do

         write(6,*)'pnucleon ', pne(1,2), pne(2,2), pne(3,2), 
     +            sqrt(pne(1,2)**2+pne(2,2)**2+pne(3,2)**2)

C         write (6,*) ' '

      end do

      end

