C Example of how to use neutxs.F

      program getxs

      implicit none

      integer nuid,i,mode
      integer z1, a1, z2, a2
      real    e0, e, xs1, xs2

20871.6

      nuid = 14
      e0   = 20.8716
      z1    =  29
      a1    = 63
      z2    = 30
      a2    = 65  
      mode  = 0

      write (6,*) 'calling neutxs '

      do i = 1, 200
         e = e0 + (i-1.)*0.050
         call neutxs(nuid,e,z1,a1,mode,xs1)
         call neutxs(nuid,e,z2,a2,mode,xs2)
*         write (6,*) ' Cross Section: ',e,xs
         print *, e, xs1, xs2
      end do

      end

c----------------------------------------------------
