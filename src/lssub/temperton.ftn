#     if !defined (nombre_max_de_facteurs)
#         define   nombre_max_de_facteurs 20
#     endif
#     if !defined (lot_maximum)
#         define   lot_maximum 16
#     endif
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
C     $Log: temperton.ftn,v $
C     Revision 3.6  2014/09/25 18:42:04  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.5  2010/03/27 21:11:47  dugas
C     - Les routines SETFFT_RD et SET99_RD remplacent SET772.
C     - La routine FFT_RD (appellant FFT991_m8 qu'on retrouve
C       dans librmn_011.a) remplace FFT772.
C
C     Revision 3.4  2008/04/28 21:38:53  dugas
C     Corriger l'usage des macros pour r.gppf (passe 2).
C
C     Revision 3.3  2003/09/15 16:19:46  dugas
C     Re-definir le macro lot_maximum de 1024 a 4 suite a
C      une modification equivalente dans les versions libpriv.a
C      de RPASSM8 et QPASSM8 sur nos frontaux survenues recamment.
C
C     Revision 3.2  1996/09/19 11:00:31  armnrbd
C     Correction temporaire pour QPASS8 et RPASS8.
C
C     Revision 3.1  1994/11/17  14:14:10  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:56:14  13:56:14  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:32:16  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.2  93/02/16  11:26:38  armnrbd
C     Enlever qpassm et rpassm.
C     Renommer set77 et fft771 en set772 et fft772.
C     
C     Revision 1.1  92/06/17  15:54:31  armnrbd
C     Desactiver le parametre WORK dans l'appel a FFT771.
C     Allouer un champs de travail correspondant avec HPALLOC.
C     
C     Revision 1.0  92/02/21  11:34:44  armnrbd
C     Initial revision
C     
C Copyright 1981-2007 ECMWF
C 
C Licensed under the GNU Lesser General Public License which
C incorporates the terms and conditions of version 3 of the GNU
C General Public License.
C See LICENSE and gpl-3.0.txt for details.
C
C  this file contains a few modifications to the ECMWF original code 
C  introduced implicit none
C  appended _RD to the original names to avoid accidents should the new vode
C    be used together with the original code
C  introduced a few top level entry points( setfft8, ffft8, etc ...) that needed
C     no work array and/or no explicit factor/trig constants initialization
C  transforms done with 8 byte reals while input/output declared as real
C  not appropriate for OMP parrallel execution (static allocations)
C  return error codes
C
      subroutine setfft_RD( n,ierr )
      implicit none
      integer  n,ierr
      external set99_RD

      integer npts
      real *8, pointer, dimension(:) :: trigs
      integer, parameter :: maxfac=nombre_max_de_facteurs
      integer, dimension(maxfac) :: ifac
      common /QQQRD_FFFT8_DRQQQ/ trigs,ifac,npts

      data npts /-1/

      ierr=0
      if(n .eq. npts) return
      if(n .gt. npts) then
         if(npts .gt. 0) deallocate(trigs)
         allocate(trigs(n+2),stat=ierr)
      endif
      if (ierr /= 0) then
         IERR=2
         WRITE(6,'("UNABLE TO ALLOCATE WORKING TRIG ARRAY.")')
         RETURN
      endif
      npts = n
      ifac=0
      call set99_RD(trigs,ifac,npts,ierr)

      return
      end

      SUBROUTINE SET99_RD(TRIGS,IFAX,N,ierr)
      implicit none
      integer N, IFAX(N),ierr
      real  *8 TRIGS(N)
      INTEGER JFAX(10),LFAX(7)
      integer ixxx, nil,nhl,k,nu,ifac,l,nfax,i
      real  *8 del, angle

C     SUBROUTINE 'SET99' - COMPUTES FACTORS OF N & TRIGONOMETRIC
C     FUNCTIONS REQUIRED BY FFT99 & FFT991

      DATA LFAX/6,8,5,4,3,2,1/
      IXXX=1

      DEL=4.0*ASIN(1.0D0)/FLOAT(N)
      NIL=0
      NHL=(N/2)-1
      DO 10 K=NIL,NHL
      ANGLE=FLOAT(K)*DEL
      TRIGS(2*K+1)=COS(ANGLE)
      TRIGS(2*K+2)=SIN(ANGLE)
   10 CONTINUE

C     FIND FACTORS OF N (8,6,5,4,3,2; ONLY ONE 8 ALLOWED)
C     LOOK FOR SIXES FIRST, STORE FACTORS IN DESCENDING ORDER

      NU=N
      IFAC=6
      K=0
      L=1
   20 CONTINUE
      IF (MOD(NU,IFAC).NE.0) GO TO 30
      K=K+1
      JFAX(K)=IFAC
      IF (IFAC.NE.8) GO TO 25
      IF (K.EQ.1) GO TO 25
      JFAX(1)=8
      JFAX(K)=6
   25 CONTINUE
      NU=NU/IFAC
      IF (NU.EQ.1) GO TO 50
      IF (IFAC.NE.8) GO TO 20
   30 CONTINUE
      L=L+1
      IFAC=LFAX(L)
      IF (IFAC.GT.1) GO TO 20

      ierr=3
      WRITE(6,40) N
   40 FORMAT('1N =',I4,' - CONTAINS ILLEGAL FACTORS')
      RETURN

C     NOW REVERSE ORDER OF FACTORS

   50 CONTINUE
      NFAX=K
      IFAX(1)=NFAX
      DO 60 I=1,NFAX
      IFAX(NFAX+2-I)=JFAX(I)
   60 CONTINUE
      IFAX(10)=N

      WRITE(6,70) IFAX(1),(IFAX(I),I=2,IFAX(1)+1)
   70 FORMAT('   FFT_M8  NFAX=',I4,', factors ...',10I4)

      RETURN
      END

      subroutine fft_RD( a, inc, jump, n, lot, isign,ierr )
      implicit none
      integer inc, jump, n, lot, isign, ierr
      real    a(*)

      integer npts
      real *8, pointer, dimension(:) :: trigs
      integer, parameter :: maxfac=nombre_max_de_facteurs
      integer, dimension(maxfac) :: ifac
      common /QQQRD_FFFT8_DRQQQ/ trigs,ifac,npts

      external fft991_m8
C
C     CHANGE maxlot ON VECTOR MACHINES
C
      integer, parameter :: maxlot=lot_maximum
      real *8, save, allocatable, dimension(:)   :: aa
      real *8, save, allocatable, dimension(:,:) :: work
      integer, save :: nwork=-1, naa=-1
      integer i,slice

      ierr=0

      if (n /= npts) then
         call setfft_RD( n,ierr )
         if (ierr /= 0) return
      endif

      if (npts > nwork) then
         if (nwork > 0) deallocate( work )
         allocate( work(npts+2,maxlot), stat=ierr )
         if (ierr /= 0) then
            ierr = 2
            return
         endif
         nwork = npts
      endif

      if (LOT*JUMP > naa) then
         if (naa > 0) deallocate( aa )
         allocate( aa(LOT*JUMP), stat=ierr )
         if (ierr /= 0) then
            ierr = 2
            return
         endif
         naa = LOT*JUMP
      endif

      aa(1:LOT*JUMP) = a(1:LOT*JUMP)

      do i=1,lot,maxlot
         slice=min(maxlot,lot+1-i)
         call fft991_m8( aa(1+(i-1)*jump), work,
     %                   trigs, ifac, inc, jump, n, slice, isign)
      enddo

      a(1:LOT*JUMP) = aa(1:LOT*JUMP)

      return
      end
