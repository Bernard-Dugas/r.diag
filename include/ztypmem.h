      real*8           rrot(3,3)
      common /D_RMAT/  rrot

      real, dimension (:), pointer :: alon,alat,lonr,latr
      real             gnplon,gnplat,longpol,rlonoff
      integer          zip1,zip2,zip3

      common /ztypmem/ alon,alat,lonr,latr
      common /ztypmem/ gnplon,gnplat,longpol,rlonoff
      common /ztypmem/ zip1,zip2,zip3
