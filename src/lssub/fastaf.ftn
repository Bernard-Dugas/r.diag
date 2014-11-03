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
C     $Log: fastaf.ftn,v $
C     Revision 3.6  2014/09/25 18:42:02  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.5  1999/11/03 21:22:43  armnrbd
C     Bug mineur...
C
C     Revision 3.4  1999/11/03 20:44:39  armnrbd
C     Ajouter la routine FASTX8.
C     Accumuler en Real*8 dans STAF.
C
C     Revision 3.3  1997/05/07 18:35:50  armnrbd
C     Interchanger les boucles en ILEV et NL.
C
C     Revision 3.2  1995/11/01  20:04:03  armnrbd
C     De-allouer les appels multi-taches.
C
C     Revision 3.1  1994/11/17  14:13:10  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:55:22  13:55:22  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.1  93/10/14  21:04:38  armnrbd
C     Petite correction pour isoler le mode MP de SGI.
C     
C     Revision 2.0  93/10/13  13:31:37  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.1  93/08/19  16:20:20  armnrbd
C     Modifications cosmetiques.
C     
C     Revision 1.0  92/02/21  11:32:20  armnrbd
C     Initial revision
C     

      SUBROUTINE FAST (S,F,LSR,LM,LA,ILH,ILEV,WLP) 

C     * JUL 21/83 - R.LAPRISE. 

C     * MULTIPLE LEGENDRE TRANSFORMS FROM FOURIER TO SPECTRAL. 
C     * ACCUMULATE THE CONTRIBUTIONS OF FOURIER COEFFICIENTS 
C     * ON ONE LATITUDE CIRCLE INTO THE LATITUDE INTEGRAL 
C     * TO OBTAIN SPHERICAL HARMONIC COEFFICIENTS. 

C     * S    = ACCUMULATED SPHERICAL HARMONIC COEFFICIENTS, 
C     * F    = FOURIER COEFFICIENTS AT ONE LATITUDE, 
C     * WLP  = LEGENDRE POLYNOMIALS NORMALIZED BY GAUSSIAN WEIGHT, 
C     * LM   = NUMBER OF E-W WAVE NUMBERS (M=1,LM), 
C     * ILH  = FIRST DIMENSION OF COMPLEX F, 
C     * LA   = FIRST DIMENSION OF COMPLEX S, 
C     * ILEV = NUMBER OF LEVELS FOR WHICH TRANSFORM IS PERFORMED. 

      IMPLICIT none

      INTEGER  L,M,N, LM,LA,ILH,ILEV,LSR(2,1), NA,NB,NL,NR
      REAL     S(2,LA,ILEV),F(2,ILH,ILEV), F1,F2
      REAL*8   WLP(1), WP

CC$   INTEGER  mpserv, Bidon
CC$   EXTERNAL mpserv

C----------------------------------------------------------------------- 
CC$doacross local( M,NA,NL,NR,N,WP,L,NB,F1,F2 )

      DO 500 M=1,LM 
      NA=LSR(2,M) 
      NL=LSR(1,M) 
      NR=LSR(1,M+1)-1 

cc    IF (ILEV.GT.NR-NL+1) THEN 
cc 
cc        DO  200 N=NL,NR 
cc            WP=WLP(NA) 
cc            DO  100 L=1,ILEV 
cc                S(1,N,L)=S(1,N,L)+WP*F(1,M,L) 
cc                S(2,N,L)=S(2,N,L)+WP*F(2,M,L) 
cc100         CONTINUE 
cc            NA=NA+1 
cc200     CONTINUE 
cc
cc    ELSE 
cc
          DO  L=1,ILEV 
              NB=NA 
              F1=F(1,M,L) 
              F2=F(2,M,L) 
              DO  N=NL,NR 
                  S(1,N,L)=S(1,N,L) +WLP(NB)*F1 
                  S(2,N,L)=S(2,N,L) +WLP(NB)*F2 
                  NB=NB+1 
              END DO
          END DO

cc    END IF 

 500  CONTINUE 

CC$    Bidon = mpserv('BLOCK',Bidon)

      RETURN 
C----------------------------------------------------------------------- 

      END 
      SUBROUTINE FASTX8 (S,F,LSR,LM,LA,ILH,ILEV,WLP) 

C     * JUL 21/83 - R.LAPRISE. 

C     * MULTIPLE LEGENDRE TRANSFORMS FROM FOURIER TO SPECTRAL. 
C     * ACCUMULATE THE CONTRIBUTIONS OF FOURIER COEFFICIENTS 
C     * ON ONE LATITUDE CIRCLE INTO THE LATITUDE INTEGRAL 
C     * TO OBTAIN SPHERICAL HARMONIC COEFFICIENTS. 

C     * S    = ACCUMULATED SPHERICAL HARMONIC COEFFICIENTS, 
C     * F    = FOURIER COEFFICIENTS AT ONE LATITUDE, 
C     * WLP  = LEGENDRE POLYNOMIALS NORMALIZED BY GAUSSIAN WEIGHT, 
C     * LM   = NUMBER OF E-W WAVE NUMBERS (M=1,LM), 
C     * ILH  = FIRST DIMENSION OF COMPLEX F, 
C     * LA   = FIRST DIMENSION OF COMPLEX S, 
C     * ILEV = NUMBER OF LEVELS FOR WHICH TRANSFORM IS PERFORMED. 

      IMPLICIT none

      INTEGER  L,M,N, LM,LA,ILH,ILEV,LSR(2,1), NA,NB,NL,NR
      REAL     F(2,ILH,ILEV), F1,F2
      REAL*8   S(2,LA,ILEV),WLP(1), WP

CC$   INTEGER  mpserv, Bidon
CC$   EXTERNAL mpserv

C----------------------------------------------------------------------- 
CC$doacross local( M,NA,NL,NR,N,WP,L,NB,F1,F2 )

      DO  L=1,ILEV 

          DO 500 M=1,LM 

              NA=LSR(2,M) 
              NL=LSR(1,M) 
              NR=LSR(1,M+1)-1 

cc        IF (ILEV.GT.NR-NL+1) THEN 
cc 
cc            DO  200 N=NL,NR 
cc                WP=WLP(NA) 
cc                DO  100 L=1,ILEV 
cc                    S(1,N,L)=S(1,N,L)+WP*F(1,M,L) 
cc                    S(2,N,L)=S(2,N,L)+WP*F(2,M,L) 
cc100             CONTINUE 
cc                NA=NA+1 
cc200         CONTINUE 
cc
cc        ELSE 
cc
              NB=NA 

              F1 = F(1,M,L) 
              F2 = F(2,M,L) 

              DO  N=NL,NR 
                  S(1,N,L) = S(1,N,L)+WLP(NB)*F1
                  S(2,N,L) = S(2,N,L)+WLP(NB)*F2
                  NB       = NB+1 
              END DO
                  
cc        END IF 

 500      CONTINUE 

      END DO

CC$    Bidon = mpserv('BLOCK',Bidon)

      RETURN 
C----------------------------------------------------------------------- 

      END 
      SUBROUTINE STAF (F,S,LSR,LM,LA,ILH,ILEV,ALP) 

C     * JUL 19/83 - R.LAPRISE. 

C     * DO A LEGENDRE TRANSFORM ON SPERICAL HARMONIC COEFFICIENTS (S) 
C     * TO GET FOURIER COEFFICIENTS (F) FOR ALL LEVELS AT THIS 
C     * LATITUDE. 

C     * ALP      = VALUE OF THE LEGENDRE POLYNOMIALS AT THIS LATITUDE. 
C     * LSR(1,M) = START OF ROW M OF S. 
C     * LSR(2,M) = START OF ROW M OF ALP. 
C     *            LSR IS DIMENSIONNED (2,LM+1). 
C     * LM       = NUMBER OF ROWS IN ARRAYS S,F AND ALP,  (M=1,LM) 
C     * LA       = FIRST DIMENSION OF COMPLEX S. 
C     * ILH      = FIRST DIMENSION OF COMPLEX F. 

      IMPLICIT none

      INTEGER  L,M,N,NA,NB,NL,NR, LM,LA,ILH,ILEV,LSR(2,1) 
      REAL     S(2,LA,ILEV),F(2,ILH,ILEV)
      REAL*8   ALP(1)
      REAL*8   AC1,AC2

CC$    INTEGER  mpserv, Bidon
CC$    EXTERNAL mpserv

C----------------------------------------------------------------------- 
C     * MXMA VERSION OF THE LM MATRIX MULTIPLIES.
C
C     DO 300 M=1,LM 
C     NA = LSR(2,M) 
C     NL = LSR(1,M) 
C     LEN= LSR(1,M+1)-NL 
C     CALL MXMA( S(1,NL,1), 2*LA , 2, 
C    1           ALP(NA)  , 1    , 1, 
C    2           F(1,M,1) , 2*ILH, 1, 
C    3           ILEV     , LEN  , 1) 
C     CALL MXMA( S(2,NL,1), 2*LA , 2, 
C    1           ALP(NA)  , 1    , 1, 
C    2           F(2,M,1) , 2*ILH, 1, 
C    3           ILEV     , LEN  , 1) 
C 300 CONTINUE 
C     RETURN 
C----------------------------------------------------------------------- 
C     * FORTRAN CODE TO DO THE LM MATRIX MULTIPLIES. 
  
CC$doacross local( M,NA,NL,NR,N,L )

      DO  300 M=1,LM 
          NA = LSR(2,M) 
          NL = LSR(1,M) 
          NR = LSR(1,M+1)-1 
 
          DO L=1,ILEV 

              AC1 = 0.0
              AC2 = 0.0

              NB  = NA
              DO  N=NL,NR 
                  AC1 = AC1+ALP(NB)* S(1,N,L)
                  AC2 = AC2+ALP(NB)* S(2,N,L)
                  NB  = NB+1 
              END DO

              F(1,M,L) = AC1
              F(2,M,L) = AC2
  
          END DO

  300 CONTINUE 

CC$    Bidon = mpserv('BLOCK',Bidon)
 
      RETURN 
C----------------------------------------------------------------------- 
      END 
