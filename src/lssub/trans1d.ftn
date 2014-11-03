#     if defined (RDIAG_LICENCE)
!---------------------------------- LICENCE BEGIN -------------------------------
! R.DIAG - Diagnostic tool kit for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This code is free software; you can redistribute it and/or modify it 
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------
#     endif
C     $Log: trans1d.ftn,v $
C     Revision 3.5  2014/10/16 12:00:45  dugas
C     Modifications commandees par la compilation avec GFORTRAN et IFORT.
C
C     Revision 3.4  2014/09/25 18:42:04  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.3  2011/11/22 19:15:36  dugas
C     Revenir a SETFFT8 et FFFT8.
C
C     Revision 3.2  2007/08/29 20:06:28  dugas
C     - Modifications du-03-2007 (BDenis) : Correction for the DFT case.
C     - Utiliser ALLOCATE plutot que STKMEMW.
C
C     Revision 3.1  2003/10/24 21:05:48  dugas
C     Implementer du code compatible RS6000
C
C     Revision 3.0  2000/07/24 20:39:14  armnrbd
C     Version initiale (de Jean Cote).
C
      SUBROUTINE TRANS1D(AR,R1,R2,R3,MAXSIZE,N,NLOT,ICAS,IAXE,IWAY)

      IMPLICIT NONE

!!!   ! AUTHOR B.DENIS (with help of J.Cote and B.Dugas)

!!!   ! Modifications: 27-03-2007 : Correction for the DFT case

!!!   ! THIS ROUTINE DOES THE SPECTRAL TRANSFORMS (FORWARD AND INVERSE)
!!!   ! IT WORKS ON MULTIPLE ROW OR COLUMN AT A TIME.

!!!   ! IT CAN DO TWO KINDS OF TRANSFORM:

!!!               1 => DCT (COS SHIFTED)
!!!               2 => DFT (REGULAR SIN & COS).

      integer maxsize,n,nlot,icas,iaxe,iway
      real*8 ar(maxsize), r1(maxsize),r2(maxsize),r3(maxsize)

      integer nx,ny, nn, lot, nic

      real*8 cnorm, fac

      real*8 del, angle

      real*8 zero, half, one, two
      parameter( zero = 0.0 )
      parameter( half = 0.5 )
      parameter( one  = 1.0 )
      parameter( two  = 2.0 )

      logical fast_ft
      logical lprint,  ldct, ldft
      integer i, j, k, inc, jump
      integer incr, jumpr, nfact, lotr
      integer ndr

      ndr(i,j) =  1 + ( j - 1 ) * jumpr + i * incr

      lprint  = .false.

      ldct   = icas .eq. 4
      ldft    = icas .eq. 5

      if (lprint) then 
      print *,'*********************'
      if ( ldct ) print *,'*  qcft8, qcfft8    *'
      if ( ldft  ) print *,'*   rft8,  ffft8    *'
      print *,'*********************'
      endif


!     * select the fast or the slow fourier transform 
!     * based on factorization.
      
      nfact=n
      call ngfft( nfact )


      if(nfact.eq.n) then 
         fast_ft=.true.
         if ( ldct ) call setscqr( n  , 'QCOS' )
         if ( ldft  ) call setscqr( n  , 'REAL' )
         
         if (iaxe.eq.0) then
            write (6,*) ' '
            write (6,*)'    using fast transform in x'
         else
            write (6,*)'    using fast transform in y'
         endif
         
      else
         fast_ft=.false.
         if (iaxe.eq.0) then
            write (6,*) ' '
            write (6,*) '    using slow transform in x !!! '
         else
            write (6,*) ' '
            write (6,*) '    using slow transform in y !!! '
         endif
      endif

      cnorm = sqrt( two/n )
      del   = acos( - one )/n

      nic  = n + 1

         if     ( IAXE .eq. 0 ) then

!!!   * on veut par exemple:

!!!   appel de trans1d    liste d'arguments  appel pour
!!!   dans main           de trans1d         routine fft
!!!   =================   =================  ===========
!!!
!!!   NX=NID    -->        N      -->        n=JUMP=N=NID
!!!   NY=NJD    -->        NLOT   -->        lot=NLOT=NJD 
!!!                                          inc=1
            lot  = NLOT

            lotr  = NLOT + 2    

            inc   = 1
            jump  = n
            incr  = 1
            jumpr = n + 2
         else

!!!   appel de trans1d    liste d'arguments  appel pour
!!!   dans main           de trans1d         routine fft
!!!   =================   =================  ===========
!!!
!!!   nx=NJD    -->        N      -->           n=NJD
!!!   ny=NID    -->        NLOT   -->         lot=INC=NLOT=NID 
!!!                                          jump=1
            lot   = NLOT
            lotr  = NLOT + 2    
            inc   = NLOT  
            jump  = 1
            incr  = lotr 
            jumpr = 1
         endif

!------------------------------------------------------------------
!
!    case of qcft8, qcfft8
!
!------------------------------------------------------------------
         if ( ldct ) then

           if (iway.eq.-1) then             !  gridpoint -> spectral
               if (fast_ft) then
                  call qcfft8(ar, inc,jump,lot,-1)
               else
                  r1=0.0
                  call qcft8 (r1,ar,inc,jump,lot,-1,n )
                  ar=r1
               endif

            else                             ! spectral -> gridpoint

               if (fast_ft) then
                  call qcfft8(ar, inc,jump,lot,+1)
               else
                  r1=0.0
                  call qcft8 (r1,ar,inc,jump,lot,+1,n )
                  ar=r1
               endif
            endif

         endif
!------------------------------------------------------------------
!
! CASE OF RFT8, FFFT8
!
!------------------------------------------------------------------
         if ( ldft ) then

            if ( lprint ) then
               print *,' on imprime ar avant transforme'
               do j=1,lot
                  print *,'j = ',j,(ar(ndr(i,j)),i=0,2*(n/2)+1)
               enddo
            endif
         
            if (iway.eq.-1) then

               if (fast_ft) then                   ! gridpoint -> spectral
                  call ffft8( ar, incr,jumpr,lotr,-1) 
               else
                  r1=0.0
                  call rft8 ( r1,ar,incr,jumpr,lotr,-1,n)
                  ar=r1
               endif
               
            else                                    ! spectral -> gridpoint
               if (fast_ft) then
                  call ffft8( ar, incr,jumpr,lotr,+1) 
               else
                  r1=0.0
                  call rft8 ( r1,ar,incr,jumpr,lotr,+1,n)
                  ar=r1
               endif
            endif
            
            if ( lprint ) then
               print *,' On imprime ar apres transforme'
               do j=1,lot
                  print *,'j = ',j,(ar(ndr(i,j)),i=0,2*(n/2)+1)
               enddo
            endif

         endif

      RETURN 
      END

C     subroutine 'setscqr' - sets up common 'comfft8x' required by
C                            cfft8, sfft8, qcfft8, qsfft8
C                                       and
C                            calls setfft8 to set up 'comfft8'
C                            required by ffft8 which is used by 
C                            the 4 previous transforms
C
Cauthor jean cote Sept 99
C
Carguments
C   i      - nf        - number of grid points (length of transform)
C   i      - case      - coded name of desired transform
C                            case = 'COS'  for real cosine         (cfft8)
C                            case = 'SIN'  for real sine           (sfft8)
C                            case = 'QCOS' for real shifted cosine (qcfft8)
C                            case = 'QSIN' for real shifted sine   (qsfft8)
C                            case = 'REAL' for real periodic       (ffft8)
C
CC
C
      subroutine setscqr( nf, case )

      implicit none

      integer nf
      character*(*) case

      integer n, m, nstore
      REAL*8, POINTER :: SSIN(:), CCOS(:), QSIN(:)
      COMMON / COMFFT8X / SSIN, CCOS, QSIN, N, M, NSTORE
C
      integer    i, ier
      real *8    del, angle
      character  alloue*17
C
      character(4) icase
      real*8 zero, half, one, two, four
      parameter( zero  = 0.0 )
      parameter( half  = 0.5 )
      parameter( one   = 1.0 )
      parameter( two   = 2.0 )
      parameter( four  = 4.0 )
C
      data alloue/'PAS ENCORE ALLOUE'/

      IF (ALLOUE.NE.'PAS ENCORE ALLOUE') DEALLOCATE(  SSIN,CCOS,QSIN )
      IF (ALLOUE.EQ.'PAS ENCORE ALLOUE') ALLOUE = 'DEJA ALLOUE'

C     n  =  length of auxiliary real periodic fourier transform (ffft8)

      icase = case
      call low2up( icase,icase )

      if     ( icase .eq. 'SIN' ) then
         n = nf + 1
      elseif ( icase .eq. 'COS' ) then
         n = nf - 1
      elseif ( icase .eq. 'REAL'   .or.
     %         icase .eq. 'QSIN'   .or.
     %         icase .eq. 'QCOS' ) then
         n = nf
      else
         print *,'ERROR in SETSCQR -> case = ',icase
         return
      endif
      
      m      = n/2
      nstore = n + 2

      ALLOCATE( SSIN(N-M-1),CCOS(N-M-1),QSIN(0:M-1) )

      del = acos( - one )/n

      do i=1,n-m-1

         angle = i * del
         ccos( i ) = cos( angle )
         ssin( i ) = sin( angle )

      enddo

      do i=0,m-1

         qsin( i ) = sin( ( i + half ) * del )

      enddo

      call setfft8( n )

      return
      end
C     subroutine 'cfft8' - multiple fast real cosine transform

C     real cosine transform of length n+1
C     the transform is its own inverse
C     self inverting implementation of Numerical Recipes pp. 508-512

C     created: Sept 99 by j.cote, rpn

C     a     is the array containing input and output data
C     inc   is the increment within each data 'vector'
C          (e.g. inc=1 for consecutively stored data)
C     jump  is the increment between the start of each data vector
C     lot   is the number of data vectors

C     definition of transform:
C     -------------------------

C     r(k) = sqrt(2/n)*sum(i=1,...,n-1)(a(i)*cos(i*k*pi/n))
C           +sqrt(1/(2n))*( a(0) + (-1)**k * a(n) )
C
C     Note for 'a' stored as a(n1,n2) then
C
C        for a transform along the first dimension
C
C           inc   = 1
C           jump  = n1
C
C        for a transform along the second dimension
C
C           inc   = n1
C           jump  = 1
C
C     The subroutine SETSCQR must have been called to set-up
C     the commons COMFFT8 and COMFFT8X
C
C-----------------------------------------------------------------------
C
      subroutine cfft8( a, inc, jump, lot )

      implicit none

      integer inc, jump, lot
      real*8  a(*)

      real*8  ai, as, ya, ys, xnor
      integer i, j, k, is

      integer n, m, nstore
      REAL*8, POINTER :: SSIN(:), CCOS(:), QSIN(:)
      COMMON / COMFFT8X / SSIN, CCOS, QSIN, N, M, NSTORE

      integer j0, jlot

      REAL*8,  DIMENSION (:), ALLOCATABLE :: W

      real*8 zero, half, one, two, four
      parameter( zero = 0.0 )
      parameter( half = 0.5 )
      parameter( one  = 1.0 )
      parameter( two  = 2.0 )
      parameter( four = 4.0 )

      integer ija, ijw
      ija(i,j) = 1 + (j0+j-1)*jump + i*inc 
      ijw(i,j) = j + i*511 
C 
C     allocate w work array
C
      ALLOCATE( W(511*NSTORE) )

      xnor = sqrt( half * n )

      do 100 j0=0,lot-1,511
      jlot  = min( 511, lot - j0 )

      i  = 1
      is = n - i
      do j=1,jlot
         ai= a( ija(i ,j) )
         as= a( ija(is,j) )
         ya = ( as - ai )
         a( ija( 1,j) ) = - ya * ccos( i )
         ya = two *  ssin( i ) * ya * xnor
         ys = ( as + ai ) * xnor 
         w( ijw(i ,j) ) = ys + ya
         w( ijw(is,j) ) = ys - ya
         ys = ( a( ija(n,j) ) + a( ija(0,j) ) ) * xnor
         a( ija(1,j) ) = a( ija(1,j) ) -
     %                   half * ( a( ija(n,j) ) - a( ija(0,j) ) )
         w( ijw(0,j) ) = ys
         w( ijw(n,j) ) = ys
      enddo

      do i=2,n-m-1

         is = n - i
         do j=1,jlot
            ai= a( ija(i ,j) )
            as= a( ija(is,j) )
            ya = ( as - ai )
            a( ija( 1,j) ) = a( ija(1,j) ) - ya * ccos( i )
            ya = two *  ssin( i ) * ya * xnor
            ys = ( as + ai ) * xnor 
            w( ijw(i ,j) ) = ys + ya
            w( ijw(is,j) ) = ys - ya
         enddo

      enddo

      do j=1,jlot
         a( ija(1,j) ) = a( ija(1,j) )/xnor
      enddo

      if ( n .eq. 2 * m ) then
         do j=1,jlot
            w( ijw( m,j) ) = two * a( ija( m,j) ) * xnor
         enddo
      endif

      call ffft8( w, 511, 1, jlot, -1 )

      do j=1,jlot
         a( ija(0,j) ) = w( ijw(0,j) )
      enddo

      do k=1,m-1
         do j=1,jlot
            a( ija(2*k  ,j) ) = w( ijw(2*k  ,j) )
            a( ija(2*k+1,j) ) = a( ija(2*k-1,j) ) - w( ijw(2*k+1,j) )
         enddo
      enddo

      do j=1,jlot
         a( ija(2*m,j) ) = w( ijw(2*m,j) )
      enddo

      if ( n .ne. 2 * m ) then
         do j=1,jlot
            a( ija(n,j) ) = a( ija(n-2,j) ) - w( ijw(n,j) )
         enddo
      endif

  100 continue

C     deallocate w work array

      DEALLOCATE( W )

      return
      end
C     subroutine 'sfft8' - multiple fast real sine transform

C     real sine transform of length n-1
C     the transform is its own inverse
C     self inverting implementation of Numerical Recipes pp. 508-512

C     created: Sept 99 by j.cote, rpn

C     a     is the array containing input and output data
C     inc   is the increment within each data 'vector'
C          (e.g. inc=1 for consecutively stored data)
C     jump  is the increment between the start of each data vector
C     lot   is the number of data vectors

C     definition of transform:
C     -------------------------

C     r(k) = sqrt(2/n)*sum(i=1,...,n-1)(a(i)*sin(i*k*pi/n))
C
C     Note for 'a' stored as a(n1,n2) then
C
C        for a transform along the first dimension
C
C           inc   = 1
C           jump  = n1
C
C        for a transform along the second dimension
C
C           inc   = n1
C           jump  = 1
C
C     The subroutine SETSCQR must have been called to set-up
C     the commons COMFFT8 and COMFFT8X
C
C-----------------------------------------------------------------------
C
      subroutine sfft8( a, inc, jump, lot )

      implicit none

      integer inc, jump, lot
      real*8  a(*)

      real*8  ai, as, ya, ys, xnor
      integer i, j, k, is

      integer n, m, nstore
      REAL*8, POINTER :: SSIN(:), CCOS(:), QSIN(:)
      COMMON / COMFFT8X / SSIN, CCOS, QSIN, N, M, NSTORE

      integer j0, jlot

      REAL*8,  DIMENSION (:), ALLOCATABLE :: W

      real*8 zero, half, one, two, four
      parameter( zero = 0.0 )
      parameter( half = 0.5 )
      parameter( one  = 1.0 )
      parameter( two  = 2.0 )
      parameter( four = 4.0 )

      integer ija, ijw
      ija(i,j) = 1 + (j0+j-1)*jump + (i-1)*inc 
      ijw(i,j) = j + i*511 

C     allocate w work array

      ALLOCATE( W(511*NSTORE) )

      xnor = sqrt( half * n )

      do 100 j0=0,lot-1,511
      jlot  = min( 511, lot - j0 )

      do j=1,jlot
         w( ijw(0,j) ) = zero
      enddo

      do i=1,n-m-1

         is = n - i
         do j=1,jlot
            ai= a( ija(i ,j) )
            as= a( ija(is,j) )
            ya=  ( as - ai ) * xnor
            ys=  two * ssin( i ) * ( as + ai ) * xnor
            w( ijw(i ,j) ) = ys + ya
            w( ijw(is,j) ) = ys - ya
         enddo

      enddo

      if ( n .eq. 2 * m ) then
         do j=1,jlot
            w( ijw(m,j) ) = four * a( ija(m,j) ) * xnor
         enddo
      endif
 
      call ffft8( w, 511, 1, jlot, -1 )
 
      do j=1,jlot
         a( ija(1,j) ) = half * w( ijw(0,j) )
      enddo

      do k = 2 , n-2 , 2

         do j=1,jlot
            a( ija(k+1,j) ) = w( ijw(k,j)   ) + a( ija(k-1,j) )
            a( ija(k,  j) ) = w( ijw(k+1,j) )
         enddo
      enddo

      if ( n .ne. 2 * m ) then
         do j=1,jlot
            a( ija(n-1,j) ) = w( ijw(n,j) )
         enddo
      endif

  100 continue

C     deallocate w work array

      DEALLOCATE( W )

      return
      end
C     subroutine 'qcfft8' - multiple fast real shifted cosine transform

C     real shifted cosine transform of length n
C     implementation inspired by Numerical Recipes pp. 513
C     but with a different normalization

C     created: Sept 99 by j.cote, rpn

C     a     is the array containing input and output data
C     inc   is the increment within each data 'vector'
C          (e.g. inc=1 for consecutively stored data)
C     jump  is the increment between the start of each data vector
C     lot   is the number of data vectors
C     isign = +1 for transform from fourier to gridpoint
C           = -1 for transform from gridpoint to fourier

C     definition of transform:
C     -------------------------

C     isign=+1: r(i) = sum(k=0,...,n-1)(a(k)*cos((i+1/2)*k*pi/n))

C     isign=-1: r(k) = sum(i=0,...,n-1)(a(i)*cos((i+1/2)*k*pi/n))
C                      * ((2-delta(k,0))/n)

C     Note for 'a' stored as a(n1,n2) then
C
C        for a transform along the first dimension
C
C           inc   = 1
C           jump  = n1
C
C        for a transform along the second dimension
C
C           inc   = n1
C           jump  = 1
C
C     The subroutine SETSCQR must have been called to set-up
C     the commons COMFFT8 and COMFFT8X
C
C-----------------------------------------------------------------------
C
      subroutine qcfft8( a, inc, jump, lot, isign )

      implicit none

      integer inc, jump, lot, isign
      real*8  a(*)

      real*8  ai, as, ya, ys, c, s, rr
      integer i, j, k, is, k1, kk

      integer n, m, nstore
      REAL*8, POINTER :: SSIN(:), CCOS(:), QSIN(:)
      COMMON / COMFFT8X / SSIN, CCOS, QSIN, N, M, NSTORE

      integer j0, jlot

      REAL*8,  DIMENSION (:), ALLOCATABLE :: W

      real*8 zero, half, one, two, four
      parameter( zero = 0.0 )
      parameter( half = 0.5 )
      parameter( one  = 1.0 )
      parameter( two  = 2.0 )
      parameter( four = 4.0 )

      integer ija, ijw
      ija(i,j) = 1 + (j0+j-1)*jump + i*inc
      ijw(i,j) = j + i*511

C     allocate w work array

      ALLOCATE( W(511*NSTORE) )

      do 100 j0=0,lot-1,511
      jlot  = min( 511, lot - j0 )

      if ( isign .eq. -1 ) then

C     transform from gridpoint to Fourier

         do i = 0 , m-1

            is = n - i - 1

            do j=1,jlot
               ai = a( ija(i ,j) )
               as = a( ija(is,j) )
               ys = ai + as
               ya =  two *  qsin( i ) * ( as - ai )
               w( ijw(i ,j) ) = ys + ya
               w( ijw(is,j) ) = ys - ya
            enddo

         enddo

         if ( n .ne. 2 * m )  then
            do j=1,jlot
               w( ijw(m,j) ) = two * a( ija(m,j) )
            enddo
         endif

         call ffft8( w, 511, 1, jlot, -1 )

         do j=1,jlot
            a( ija(0,j) ) = w( ijw(0,j) ) * half
         enddo

         do k = 1 , m

            kk = 2*k
            k1 = kk + 1

            if ( k .lt. m .or. n .ne. 2 * m ) then
               c = ccos( k )
               s = ssin( k )
               do j=1,jlot
                  a( ija(kk-1,j) ) = -s*w( ijw(kk,j) )+c*w( ijw(k1,j) )
                  a( ija(kk  ,j) ) =  c*w( ijw(kk,j) )+s*w( ijw(k1,j) )
               enddo
            else
               do j=1,jlot
                  a( ija(kk-1,j) ) = -w( ijw(kk,j) )
               enddo
            endif

         enddo
         if ( n .eq. 2 * m )  then
            do j=1,jlot
               a( ija(n-1,j) ) = a( ija(n-1,j) ) * half
            enddo
         endif
         do k = 2*m-3 , 1 , -2
            do j=1,jlot
               a( ija(k,j) ) = a( ija(k,j) ) + a( ija(k+2,j) )
            enddo
         enddo


      elseif ( isign .eq. +1 ) then

C     transform from Fourier to gridpoint

         do j=1,jlot
            w( ijw(0,j) ) = a( ija(0,j) )
            w( ijw(1,j) ) = zero
         enddo

         do k = 2 , n-1 , 2
            do j=1,jlot
               w( ijw(k,j) ) = a( ija(k,j) ) * half
            enddo
         enddo
         do k = 3 , 2*m-1 , 2
            do j=1,jlot
              w( ijw(k,j) ) = ( a( ija(k-2,j) ) - a( ija(k,j) ) ) * half
            enddo
         enddo
         if ( n .eq. 2 * m )  then
            c = one
         else
            c = half
         endif
         do j=1,jlot
            w( ijw(2*m+1,j) ) = a( ija(2*m-1,j) ) * c
         enddo

         do k = 1 , m

            if ( k .lt. m .or. n .ne. 2 * m ) then
               c = ccos( k )
               s = ssin( k )
            else
               c = zero
               s = one
            endif

            kk = 2*k
            k1 = kk + 1
            do j=1,jlot
               rr = w( ijw(kk,j) )
               w( ijw(kk,j) ) =  c * rr - s * w( ijw(k1,j) )
               w( ijw(k1,j) ) =  s * rr + c * w( ijw(k1,j) )
            enddo

         enddo

         call ffft8( w, 511, 1, jlot, +1 )

         do i = 0 , m-1

            is = n - i - 1

            do j=1,jlot
               ys = ( w( ijw(i ,j) ) + w( ijw(is,j) ) ) * half
               ya = ( w( ijw(is,j) ) - w( ijw(i ,j) ) ) /
     %              ( four * qsin( i ) )
               a( ija(i ,j) ) = ys + ya
               a( ija(is,j) ) = ys - ya

            enddo

         enddo
         if ( n .ne. 2 * m )  then
            do j=1,jlot
               a( ija(m,j) ) = w( ijw(m,j) )
            enddo
         endif

      endif

  100 continue

C     deallocate w work array

      DEALLOCATE( W )

      return
      end
C     subroutine 'qsfft8' - multiple fast real shifted sine transform

C     real shifted sine transform of length n
C     implementation inspired by Numerical Recipes pp. 513
C     but with a different normalization

C     created: Sept 99 by j.cote, rpn

C     a     is the array containing input and output data
C     inc   is the increment within each data 'vector'
C          (e.g. inc=1 for consecutively stored data)
C     jump  is the increment between the start of each data vector
C     lot   is the number of data vectors
C     isign = +1 for transform from fourier to gridpoint
C           = -1 for transform from gridpoint to fourier

C     definition of transform:
C     -------------------------

C     isign=+1: r(i) = sum(k=1,...,n)(a(k)*sin((i+1/2)*k*pi/n))

C     isign=-1: r(k) = sum(i=0,...,n-1)(a(i)*sin((i+1/2)*k*pi/n))
C                      * ((2-delta(k,n))/n)

C     Note for 'a' stored as a(n1,n2) then
C
C        for a transform along the first dimension
C
C           inc   = 1
C           jump  = n1
C
C        for a transform along the second dimension
C
C           inc   = n1
C           jump  = 1
C
C     The subroutine SETSCQR must have been called to set-up
C     the commons COMFFT8 and COMFFT8X
C
C-----------------------------------------------------------------------
C
      subroutine qsfft8( a, inc, jump, lot, isign )

      implicit none

      integer inc, jump, lot, isign
      real*8  a(*)

      real*8  ai, as, ya, ys, c, s, rr
      integer i, j, k, is, k1, kk

      integer n, m, nstore
      REAL*8, POINTER :: SSIN(:), CCOS(:), QSIN(:)
      COMMON / COMFFT8X / SSIN, CCOS, QSIN, N, M, NSTORE

      integer j0, jlot

      REAL*8,  DIMENSION (:), ALLOCATABLE :: W

      real*8 zero, half, one, two, four
      parameter( zero = 0.0 )
      parameter( half = 0.5 )
      parameter( one  = 1.0 )
      parameter( two  = 2.0 )
      parameter( four = 4.0 )

      integer ijg, ijf, ijw
      ijg(i,j) = 1 + (j0+j-1)*jump +     i*inc 
      ijf(i,j) = 1 + (j0+j-1)*jump + (i-1)*inc 
      ijw(i,j) = j + i*511 

C     allocate w work array

      ALLOCATE( W(511*NSTORE) )

      do 100 j0=0,lot-1,511
      jlot  = min( 511, lot - j0 )

      if ( isign .eq. -1 ) then
C
C     transform from gridpoint to Fourier
C
         do i=0,m-1

            is = n - i - 1

            do j=1,jlot
               ai = a( ijg(i ,j) )
               as = a( ijg(is,j) )
               ya = as - ai
               ys =  two *  qsin( i ) * ( as + ai )
               w( ijw(i ,j) ) = ys + ya
               w( ijw(is,j) ) = ys - ya
            enddo

         enddo
         if ( n .ne. 2 * m )  then
            do j=1,jlot
               w( ijw(m,j) ) = four * a( ijg(m,j) )
            enddo
         endif

         call ffft8( w, 511, 1, jlot, -1 )

         do k=1,m

            if ( k .lt. m .or. n .ne. 2*m ) then
               c = ccos( k )
               s = ssin( k )
            else
               c = zero
               s = one
            endif
            kk = 2 * k
            k1 = kk + 1
            do j=1,jlot
               rr = w( ijw(kk,j) )
               w( ijw(kk,j) ) =   c * rr + s * w( ijw(k1,j) )
               w( ijw(k1,j) ) = - s * rr + c * w( ijw(k1,j) )
            enddo

         enddo

         do j=1,jlot
            a( ijf(1,j) ) = half * w( ijw(0,j) )
         enddo

         do k=2,n-1,2

            do j=1,jlot
               a( ijf(k+1,j) ) = w( ijw(k  ,j) ) + a( ijf(k-1,j) )
               a( ijf(k  ,j) ) = w( ijw(k+1,j) )
            enddo

         enddo

         if ( n .eq. 2 * m ) then
            do j=1,jlot
               a( ijf(n,j) ) = w( ijw(n+1,j) )
            enddo
         endif
         do j=1,jlot
            a( ijf(n,j) ) = half * a( ijf(n,j) )
         enddo

      elseif ( isign .eq. +1 ) then
C
C     transform from Fourier to gridpoint
C
         do j=1,jlot
            w( ijw(0,j) ) = a( ijf(1,j) )
            w( ijw(1,j) ) = zero
         enddo

         do k=2,n-1,2

            do j=1,jlot
               w( ijw(k  ,j) ) = half * ( a( ijf(k+1,j) ) -
     %                                    a( ijf(k-1,j) ) )
               w( ijw(k+1,j) ) = half * a( ijf(k,j) )
            enddo

         enddo

         if ( n .eq. 2 * m ) then
            do j=1,jlot
              w( ijw(n+1,j) ) = a( ijf(n,j) )
            enddo
         else
            do j=1,jlot
              w( ijw(n-1,j) ) = w( ijw(n-1,j) ) + half * a( ijf(n,j) )
            enddo
         endif

         do k = 1 , m

            if ( k .lt. m .or. n .ne. 2*m ) then
               c = ccos( k )
               s = ssin( k )
            else
               c = zero
               s = one
            endif
            kk = 2 * k
            k1 = kk + 1

            do j=1,jlot

               rr = w( ijw(kk,j) )
               w( ijw(kk,j) ) = c * rr - s * w( ijw(k1,j) )
               w( ijw(k1,j) ) = s * rr + c * w( ijw(k1,j) )

            enddo

         enddo

         call ffft8( w, 511, 1, jlot, +1 )

         do i=0,m-1

            is = n - i - 1

            do j=1,jlot

               ys = ( w( ijw(i ,j) ) + w( ijw(is,j) ) )/
     %              ( four * qsin( i ) )
               ya = ( w( ijw(is,j) ) - w( ijw(i ,j) ) ) * half
               a( ijg(i, j) ) = ys + ya
               a( ijg(is,j) ) = ys - ya

            enddo

         enddo

         if ( n .ne. 2 * m ) then
            do j=1,jlot
               a( ijg(m,j) ) = w( ijw(m,j) ) * half
            enddo
         endif

      endif

  100 continue

C     deallocate w work array

      DEALLOCATE( W )

      return
      end
CCCs/r ngfft - calcul du prochain entier dans la suite 4, 6, 8
C              ou qui se factorise en 2, 3, et 5 quand > que 8
C              pour fft771 et ffft8, s(q)fft8, c(q)fft8
C
      subroutine ngfft( n )

      implicit none
      integer n
C
Cauteur jean cote - 1990
C
Crevision jean cote - Sept 1999, ajoute 4, 6, 8
C
Carguments
C   io     - n       - en sortie le plus petit entier >= n qui factorise
C
Cparametres
      integer l
      parameter ( l = 3 )
      integer k( l )
      data k / 2 , 3 , 5 /

      integer i,j

      if ( n .le. 8 ) then
         if ( n .le. 4 ) n = 4
         if ( n .eq. 5 ) n = 6
         if ( n .eq. 7 ) n = 8
         return
      endif

      n = n - 1
    1 n = n + 1
      i = n
    2 do 3 j=1,l
         if( mod(i,k(j)) .eq. 0 ) go to 4
    3 continue
      go to 1
    4 i = i/k(j)
      if( i .ne. 1 ) go to 2
      return
      end

C     subroutine 'cft8' - multiple slow real cosine transform

C     real cosine transform of length n+1
C     the transform is its own inverse
C     self inverting implementation of Numerical Recipes pp. 508-512

C     created: mar 10/99 by j.cote, rpn

C     r     is the array containing output data
C     a     is the array containing input data
C     inc   is the increment within each data 'vector'
C          (e.g. inc=1 for consecutively stored data)
C     jump  is the increment between the start of each data vector
C     lot   is the number of data vectors

C     definition of transform:
C     -------------------------

C     r(k) = sqrt(2/n)*sum(i=1,...,n-1)(a(i)*cos(i*k*pi/n))
C           +sqrt(1/(2n))*( a(0) + (-1)**k * a(n) )
C
C     Note for 'a' stored as a(n1,n2) then
C
C        for a transform along the first dimension
C
C           inc   = 1
C           jump  = n1
C
C        for a transform along the second dimension
C
C           inc   = n1
C           jump  = 1
C
C-----------------------------------------------------------------------
C
      subroutine cft8( r, a, inc, jump, lot, n )

C
      implicit none
C
      integer inc, jump, lot, n
      real*8  r(*), a(*)
C
      integer i, j, k
C
      real*8 zero, half, one, two
      parameter( zero = 0.0 )
      parameter( half = 0.5 )
      parameter( one  = 1.0 )
      parameter( two  = 2.0 )

      real*8 cnorm, del, angle, cc, cnorm2, faz
      integer ij
      ij(i,j) = 1 + ( j - 1 ) * jump + i * inc

      cnorm = sqrt( two/n )
      cnorm2 = sqrt( one/(2*n) )
      del   = acos( - one )/n
C
C     The transform  is along the i-th direction
C     and the indexing is 1+(j-1)*jump+i*inc (i=0,n,j=1,lot)
C     in the input & output fields
C
      faz = - one
      do k=0,n
         faz = - faz

         do j=1,lot
            r( ij(k,j) ) = cnorm2 * ( a( ij(0,j) )
     %                                  + faz * a( ij(n,j) ) )
         enddo

         do i=1,n-1

            angle = k * i * del
            cc = cnorm * cos( angle )

            do j=1,lot
               r( ij(k,j) ) = r( ij(k,j) ) + cc * a( ij(i,j) )
            enddo

         enddo

      enddo
         
      return
      end
C     subroutine 'sft8' - multiple slow real sine transform

C     real sine transform of length n-1
C     the transform is its own inverse
C     self inverting implementation of Numerical Recipes pp. 508-512

C     created: mar 10/99 by j.cote, rpn

C     r     is the array containing output data
C     a     is the array containing input data
C     inc   is the increment within each data 'vector'
C          (e.g. inc=1 for consecutively stored data)
C     jump  is the increment between the start of each data vector
C     lot   is the number of data vectors

C     definition of transform:
C     -------------------------

C     r(k) = sqrt(2/n)*sum(i=1,...,n-1)(a(i)*sin(i*k*pi/n))
C
C     Note for 'a' stored as a(n1,n2) then
C
C        for a transform along the first dimension
C
C           inc   = 1
C           jump  = n1
C
C        for a transform along the second dimension
C
C           inc   = n1
C           jump  = 1
C
C-----------------------------------------------------------------------
C
      subroutine sft8( r, a, inc, jump, lot, n )

C
      implicit none
C
      integer inc, jump, lot, n
      real*8  r(*), a(*)
C
      integer i, j, k
C
      real*8 zero, half, one, two
      parameter( zero = 0.0 )
      parameter( half = 0.5 )
      parameter( one  = 1.0 )
      parameter( two  = 2.0 )

      real*8 cnorm, del, angle, cs
      integer ij
      ij(i,j) = 1 + ( j - 1 ) * jump + ( i - 1 ) * inc

      cnorm = sqrt( two/n )
      del   = acos( - one )/n
C
C     The transform  is along the i-th direction
C     and the indexing is 1+(j-1)*jump+(i-1)*inc
C     in the input field a
C
      do k=1,n-1

         do j=1,lot
            r( ij(k,j) ) = zero
         enddo

         do i=1,n-1

            angle = k * i * del
            cs = cnorm * sin( angle )

            do j=1,lot
               r( ij(k,j) ) = r( ij(k,j) ) + cs * a( ij(i,j) )
            enddo

         enddo

      enddo
         
      return
      end
C     subroutine 'qcft8' - multiple slow real shifted cosine transform

C     real shifted cosine transform of length n
C     implementation with Clive Temperton's normalization

C     created: sept 1/99 by j.cote, rpn

C     r     is the array containing output data
C     a     is the array containing input data
C     inc   is the increment within each data 'vector'
C          (e.g. inc=1 for consecutively stored data)
C     jump  is the increment between the start of each data vector
C     lot   is the number of data vectors
C     isign = +1 for transform from fourier to gridpoint
C           = -1 for transform from gridpoint to fourier

C     definition of transform:
C     -------------------------

C     isign=+1: r(i) = sum(k=0,...,n-1)(a(k)*cos((i+1/2)*k*pi/n))

C     isign=-1: r(k) = sum(i=0,...,n-1)(a(i)*cos((i+1/2)*k*pi/n))
C                      * ((2-delta(k,0))/n)

C     Note for 'a' stored as a(n1,n2) then
C
C        for a transform along the first dimension
C
C           inc   = 1
C           jump  = n1
C
C        for a transform along the second dimension
C
C           inc   = n1
C           jump  = 1
C
C-----------------------------------------------------------------------
C
      subroutine qcft8( r, a, inc, jump, lot, isign, n )

C
      implicit none
C
      integer inc, jump, lot, isign, n
      real*8  r(*), a(*)
C
      integer i, j, k
C
      real*8 zero, half, one, two
      parameter( zero = 0.0 )
      parameter( half = 0.5 )
      parameter( one  = 1.0 )
      parameter( two  = 2.0 )

      real*8 del, angle, cc, fac
      integer ij
      ij(i,j) = 1 + ( j - 1 ) * jump + i * inc

      del   = acos( - one )/n
C
C     The transform  is along the i-th direction
C     and the indexing is 1+(j-1)*jump+ i*inc (i=0,n-1,j=1,lot)
C                      or 1+(j-1)*jump+ k*inc (k=0,n-1,j=1,lot)
C     in gridpoint and Fourier space respectively
C
      if ( isign .eq. -1 ) then
C
C     transform from gridpoint to fourier
C
         fac = one/n

         do k=0,n-1

            do j=1,lot
               r( ij(k,j) ) = zero
            enddo

            angle =  k * del

            do i=0,n-1

               cc =  cos( ( i + half ) * angle )

               do j=1,lot
                  r( ij(k,j) ) = r( ij(k,j) ) + cc * a( ij(i,j) )
               enddo

            enddo

            do j=1,lot
               r( ij(k,j) ) = fac * r( ij(k,j) )
            enddo

            fac = two/n

         enddo
         
      elseif( isign .eq. +1 ) then
C
C     transform from Fourier to gridpoint
C
         do i=0,n-1

            do j=1,lot
               r( ij(i,j) ) = a( ij(0,j) )
            enddo
 
            angle = ( i + half ) * del

            do k=1,n-1

               cc =  cos( k * angle )

               do j=1,lot
                  r( ij(i,j) ) = r( ij(i,j) ) + cc * a( ij(k,j) )
               enddo

            enddo

         enddo

      endif

      return
      end
C     subroutine 'qsft8' - multiple slow real shifted sine transform

C     real shifted sine transform of length n
C     implementation with Clive Temperton's normalization

C     created: sept 1/99 by j.cote, rpn

C     r     is the array containing output data
C     a     is the array containing input data
C     inc   is the increment within each data 'vector'
C          (e.g. inc=1 for consecutively stored data)
C     jump  is the increment between the start of each data vector
C     lot   is the number of data vectors
C     isign = +1 for transform from fourier to gridpoint
C           = -1 for transform from gridpoint to fourier

C     definition of transform:
C     -------------------------

C     isign=+1: r(i) = sum(k=1,...,n)(a(k)*sin((i+1/2)*k*pi/n))

C     isign=-1: r(k) = sum(i=0,...,n-1)(a(i)*sin((i+1/2)*k*pi/n))
C                      * ((2-delta(k,n))/n)

C     Note for 'a' stored as a(n1,n2) then
C
C        for a transform along the first dimension
C
C           inc   = 1
C           jump  = n1
C
C        for a transform along the second dimension
C
C           inc   = n1
C           jump  = 1
C
C-----------------------------------------------------------------------
C
      subroutine qsft8( r, a, inc, jump, lot, isign, n )

C
      implicit none
C
      integer inc, jump, lot, isign, n
      real*8  r(*), a(*)
C
      integer i, j, k, km
C
      real*8 zero, half, one, two
      parameter( zero = 0.0 )
      parameter( half = 0.5 )
      parameter( one  = 1.0 )
      parameter( two  = 2.0 )

      real*8 del, angle, ss, fac
      integer ij
      ij(i,j) = 1 + ( j - 1 ) * jump + i * inc

      del   = acos( - one )/n
C
C     The transform  is along the i-th direction
C     and the indexing is 1+(j-1)*jump+    i*inc (i=0,n-1,j=1,lot)
C                      or 1+(j-1)*jump+(k-1)*inc (k=1,  n,j=1,lot)
C     in gridpoint and Fourier space respectively
C
      if ( isign .eq. -1 ) then
C
C     transform from gridpoint to fourier
C
         fac = one/n

         do k=n,1,-1
            km = k - 1

            do j=1,lot
               r( ij(km,j) ) = zero
            enddo

            angle =  k * del

            do i=0,n-1

               ss =  sin( ( i + half ) * angle )

               do j=1,lot
                  r( ij(km,j) ) = r( ij(km,j) ) + ss * a( ij(i,j) )
               enddo

            enddo

            do j=1,lot
               r( ij(km,j) ) = fac * r( ij(km,j) )
            enddo

            fac = two/n

         enddo
         
      elseif( isign .eq. +1 ) then
C
C     transform from Fourier to gridpoint
C
         do i=0,n-1

            do j=1,lot
               r( ij(i,j) ) = zero
            enddo
 
            angle = ( i + half ) * del

            do k=1,n
               km = k - 1

               ss =  sin( k * angle )

               do j=1,lot
                  r( ij(i,j) ) = r( ij(i,j) ) + ss * a( ij(km,j) )
               enddo

            enddo

         enddo

      endif

      return
      end
C     subroutine 'rft8' - multiple slow real transform

C     real transform of length n

C     created: mar 12/99 by j.cote, rpn

C     r     is the array containing output data
C     a     is the array containing input data
C     inc   is the increment within each data 'vector'
C          (e.g. inc=1 for consecutively stored data)
C     jump  is the increment between the start of each data vector
C     lot   is the number of data vectors
C     isign = +1 for transform from fourier to gridpoint
C           = -1 for transform from gridpoint to fourier

C     ordering of fourier coefficients:
C         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
C         where b(0)=0, possibly b(n/2)=0 ; (2*[n/2]+2) locations required

C     ordering of gridpoint data:
C         x(0),x(1),x(2),...,x(n-1), ? , ? ; (2*[n/2]+2) locations required

C     n must be composed of factors 2,3 & 5 but does not have to be even

C     definition of transforms:
C     -------------------------

C     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
C         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
C               using the symmetry the sum restricted to k = 0, [n/2]

C     isign=-1: a(k)= (1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
C               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
C               k = 0, [n/2]

C
C     Note for 'a' stored as a(n1,n2) then
C
C        for a transform along the first dimension
C
C           inc   = 1
C           jump  = n1
C
C        for a transform along the second dimension
C
C           inc   = n1
C           jump  = 1
C
C-----------------------------------------------------------------------
C
      subroutine rft8( r, a, inc, jump, lot, isign, n )

C
      implicit none
C
      integer inc, jump, lot, isign, n
      real*8  r(*), a(*)
C
      integer i, j, k, kk, k1
C
      real*8 zero, half, one, two
      parameter( zero = 0.0 )
      parameter( half = 0.5 )
      parameter( one  = 1.0 )
      parameter( two  = 2.0 )

      real*8 oneon, del, angle, cc, ss, faz
      integer ij
      ij(i,j) = 1 + ( j - 1 ) * jump + i * inc

      oneon = one/n
      del   = acos( - one ) * two/n
C
C     The transform  is along the i-th direction
C     and the indexing is 1+(j-1)*jump+i*inc (i=0,2*[n/2]+1,j=1,lot)
C     in the input & output fields
C
      if ( isign .eq. -1 ) then
C
C     transform from gridpoint to fourier
C
      do k=0,2*(n/2)+1
         do j=1,lot
            r( ij(k,j) ) = zero
         enddo
      enddo

      do k=0,n/2

         kk = 2 * k
         k1 = kk + 1

         do i=0,n-1

            angle = k * i * del
            cc =   oneon * cos( angle )
            ss = - oneon * sin( angle )

            do j=1,lot
               r( ij(kk,j) ) = r( ij(kk,j) ) + cc * a( ij(i,j) )
               r( ij(k1,j) ) = r( ij(k1,j) ) + ss * a( ij(i,j) )
            enddo

         enddo

      enddo

      elseif ( isign .eq. +1 ) then
C
C    transform from fourier to gridpoint
C
      faz = - one * mod( n + 1, 2 )
      do i=0,n-1

         faz = - faz
         do j=1,lot
            r( ij(i,j) ) = a( ij(0,j) ) + faz * a( ij(2*(n/2),j) )
         enddo

         do k=1,(n-1)/2

            kk = 2 * k
            k1 = kk + 1
            angle = k * i * del
            cc = two * cos( angle )
            ss = two * sin( angle )

            do j=1,lot
               r( ij(i,j) ) = r( ij(i,j) ) + cc * a( ij(kk,j) )
     %                                     - ss * a( ij(k1,j) )
            enddo

         enddo

      enddo

      endif
         
      return
      end

