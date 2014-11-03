#     include "diagmacros.cdk"
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
C     $Log: fastaf2.ftn,v $
C     Revision 3.7  2014/09/25 18:31:50  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.6  2013/10/08 01:27:56  bernard
C     Ajouter STAF3 ou les variables CFc et SC sont des COMPLEX*16.
C
C     Revision 3.5  2000/03/22 19:28:25  armnrbd
C     Ajouter un include pour deinir le macro AIMAG (f77 vs f90).
C
C     Revision 3.4  2000/03/17 21:31:33  armnrbd
C     Remplacer REAL() par DBLE().
C
C     Revision 3.3  2000/03/17 03:19:14  armnrbd
C     Remplaacer IMAG par AIMAG dans FAST2 et STAF2.
C
C     Revision 3.2  1995/11/01 20:04:03  armnrbd
C     De-allouer les appels multi-taches.
C
C     Revision 3.1  1994/11/17  14:13:11  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:55:23  13:55:23  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:31:38  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.1  93/06/25  15:42:18  armnrbd
C     Implanter en mode "IMPLICIT none"
C     Declarer ALP,EPSI et WRKL en REAL*8.
C     
C     Revision 1.0  92/02/21  11:32:25  armnrbd
C     Initial revision

      SUBROUTINE fast2 (SC,LSR,LM,CFC,ALP,WEIGHT)

***    APR  1/80 - J.D.HENDERSON

***    CONTRIBUTION OF COMPLEX FOURIER COEFF IN CF TO LATITUDE 
***    INTEGRAL ADDED TO  SPHERICAL HARMONICS IN SC.  SPECTRAL 
***    HARMONICS MUST BE GLOBAL.

***    LSR(1,M) = START OF ROW M OF S
***    LSR(2,M) = START OF ROW M OF ALP. THE LONGEST ROW LENGTH 
***               PERMITTED IS 60.
***    LM       = NUMBER OF ROWS IN SC,ALP ARRAYS.
***    ALP      = ASSOCIATED LEGENDRE POLYNOMIALS FOR GIVEN LATITUDE.
***    WEIGHT   = GAUSSIAN WEIGHT FOR GIVEN LATITUDE.

      IMPLICIT   none

      COMPLEX    SC(1)
      COMPLEX*16 CFC(1)
      REAL*8     ALP(1),WEIGHT
      INTEGER    LSR(2,1),LM

      INTEGER    M,N,NA,NL,NR
      REAL       CFCR,CFCI

CC$    INTEGER    mpserv, Bidon
CC$    EXTERNAL   mpserv

C----------------------------------------------------------------------- 
CC$doacross local( M,NA,NL,NR,N,CFCR,CFCI )

      DO  220 M=1,LM 
          NA   = LSR(2,M)-1 
          NL   = LSR(1,M)
          NR   = LSR(1,M+1)-1

          CFCR = DBLE( CFC(M) )*WEIGHT
          CFCI = AIMAG( CFC(M) )*WEIGHT

          DO  220 N=NL,NR

              NA    = NA+1
              SC(N) = SC(N)+CMPLX( CFCR*ALP(NA),CFCI*ALP(NA) )

  220 CONTINUE

CC$    Bidon = mpserv('BLOCK',Bidon)

      RETURN
*---------------------------------------------------------------------

      END 
      SUBROUTINE staf2 (CFC,SC,LSR,LM,ALP)

***    APR  1/80 - J.D.HENDERSON

***    CONVERTS SPHERICAL HARMONICS IN SC TO FOURIER COEFF IN CF AT 
***    ONE PARTICULAR LATITUDE. SPECTRAL HARMONICS MUST BE GLOBAL.

***    LSR(1,M) = START OF ROW M OF S
***    LSR(2,M) = START OF ROW M OF ALP. THE LONGEST ROW LENGTH 
***               PERMITTED IS 60.
***    LM       = NUMBER OF ROWS IN SC,ALP ARRAYS.
***    ALP      = ASSOCIATED LEGENDRE POLYNOMIALS FOR GIVEN LATITUDE.

      IMPLICIT   none

      COMPLEX    SC(1)
      COMPLEX*16 CFC(1)
      REAL*8     ALP(1)
      INTEGER    LSR(2,1),LM

      REAL*8     FCR,FCI
      INTEGER    M,N,NA,NL,NR

CC$    INTEGER    mpserv, Bidon
CC$    EXTERNAL   mpserv

C----------------------------------------------------------------------- 
CC$doacross local( M,NA,NL,NR,N,FCR,FCI )

      DO  220 M=1,LM 
          NA  = LSR(2,M)-1 
          NL  = LSR(1,M)
          NR  = LSR(1,M+1)-1

          FCR = 0.
          FCI = 0.

          DO  210 N=NL,NR
              NA  = NA  + 1
              FCR = FCR + ALP(NA)*DBLE( SC(N) )
              FCI = FCI + ALP(NA)*AIMAG( SC(N) )
  210     CONTINUE

          CFC(M) = DCMPLX( FCR, FCI )

  220 CONTINUE

CC$    Bidon = mpserv('BLOCK',Bidon)

      RETURN
*----------------------------------------------------------------------

      END 
      SUBROUTINE staf3 (CFC,SC,LSR,LM)

      IMPLICIT   none

      COMPLEX*16  CFC(1),SC(1)
      INTEGER     LSR(2,1),LM

!**   APR  1/80 - J.D.HENDERSON
!     Octobre 2013 - B.Dugas : Version adaptee a Belousov pour r.diag


!**    CONVERTS SPHERICAL HARMONICS IN SC TO FOURIER COEFF IN CF AT 
!**    ONE PARTICULAR LATITUDE. SPECTRAL HARMONICS MUST BE GLOBAL.

!**    LSR(1,M) = START OF ROW M OF S
!**    LSR(2,M) = START OF ROW M OF ALP
!**    LM       = NUMBER OF ROWS IN SC,ALP ARRAYS.

      COMPLEX*16  FC
      INTEGER     M,N,NA,NL,NR

#     include  "calpi.cdk"

!----------------------------------------------------------------------- 

      DO  M=1,LM 

          NA  = LSR(1,M)-1 
          NR  = LSR(1,M+1)-LSR(1,M)

          FC = CMPLX( 0.,0. )

          DO  N=1,NR
              NA = NA + 1
              FC = FC + ALP(N,M)*SC(NA)
          END DO

          CFC(M) = FC

      END DO

      RETURN
!----------------------------------------------------------------------

      END 
