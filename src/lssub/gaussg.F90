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
!     $Log: gaussg.F90,v $
!     Revision 3.7  2015/04/21 00:35:29  dugas
!     GFORTRAN n'utilise plus DDFUN90.
!
!     Revision 3.6  2015/04/04 02:12:55  dugas
!      - Les compilateurs INTEL ne font plus appel au package DDFUN90.
!      - Renommer a gaussg.F90 pour utiliser les macros d'identification des compilateurs.
!
!     Revision 3.5  2014/09/25 18:42:02  dugas
!     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
!
!     Revision 3.4  2013/11/01 21:46:32  dugas
!     Correction a ORDLEG16 pour AIX.
!
!     Revision 3.3  2013/10/28 21:04:12  bernard
!     Ajouter GAUSSG16/ORDLEG16 et formatter pour F90.
!
!     Revision 3.2  1998/01/26 18:28:18  armnrbd
!     Remplacer ORDLEG par ORDLEG8.
!
!     Revision 3.1  1994/11/17  14:13:25  armnrbd
!     Messages informatifs quand au passage de la version 2.x a 3.1...
!     1) Les espaces en debut des noms de variables de sont plus pertinents.
!     2) Les grilles complexes de type CMPL sont maintenant supportees.
!     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
!     4) Plusieurs nouvelles cles sont disponibles au demarrage.
!
!     Revision 3.0  94/11/17  13:55:33  13:55:33  armnrbd (Bernard Dugas)
!     *** empty log message ***
!     
!     Revision 2.0  93/10/13  13:31:44  armnrbd
!     Premiere version compatible HP-UX.
!     
!     Revision 1.1  92/12/24  10:57:42  armnrbd
!     Ajouter la routine TRIGL2.
!     
!     Revision 1.0  92/02/21  11:33:03  armnrbd
!     Initial revision
!     

      SUBROUTINE gaussg16 (NZERO,F,WT,SIA,RAD,WOCS)
#     if !defined (__HOS_AIX__) && !defined (__INTEL_COMPILER_UPDATE) && !defined (__GFORTRAN__)
      USE ddmodule
#     endif
      IMPLICIT NONE
 
      INTEGER  NZERO
      real(8)  F(NZERO),WT(NZERO),SIA(NZERO),RAD(NZERO),WOCS(NZERO)

!**    THIS ROUTINE CALCULATES THE ROOTS (F) OF THE ORDINARY LEGENDRE
!**    POLYNOMIALS OF ORDER 2*NZERO.  THE FIRST STEP IS TO MAKE AN 
!**    INITIAL GUESS FOR EACH ROOT AND THEN TO USE THE ORDINARY
!**    LEGENDRE ALGORITHM (ordleg16) AND NEWTONS METHOD TO REFINE
!**    THE SOLUTION UNTIL THE CRITERION XLIM IS SATISFIED. THIS
!**    IS DONE ON THE FIRST HALF OF THE RANGE DUE TO THE SYMMETRY
!**    THE ZEROS WITH RESPECT TO F=0. THIS IS A REAL(16) VERSION
!**    IN WHICH THE EXTENDED RANGE ARITHMETIC IS PROVIDED BY THE
!**    DDFUN90 PACKAGE OF David H. Bailey (LBNL). SEE THE FOLLOWING
!**    SITE FOR MORE INFO... "crd-legacy.lbl.gov/~dhbailey/mpdist/".
!**    THE DDFUN90 VERSION IS DATED 2005-03-11.

!**       F = COSINE OF COLATITUDE 
!**      WT = CORRESPONDING GAUSSIAN WEIGHT
!**     SIA = SINE OF COLATITUDE 
!**     RAD = COLATITUDE IN RADIANS
!**    WOCS = GAUSSIAN WEIGHT / COS(COLAT)**2

      REAL(8)        ERRW8,SUMW8
      INTEGER        IR,IRP,IRM,I
#     if defined (__HOS_AIX__) || defined (__INTEL_COMPILER_UPDATE) || defined (__GFORTRAN__)
      REAL*16        FI,FI0,FI1,FN, G,GM,GP,GT, FTEMP,GTEMP, &
                     A,B,C,D, DI, DN,DN1, PI,XLIM, DOT, &
                     UN,DEUX,QUATRE,POINT5,SUMW16
#     else
      TYPE (DD_REAL) FI,FI0,FI1,FN, G,GM,GP,GT, FTEMP,GTEMP, &
                     A,B,C,D, DI, DN,DN1, PI,XLIM, DOT, &
                     UN,DEUX,QUATRE,POINT5,SUMW16
#     endif
      LOGICAL        INFO
      EXTERNAL       GET_INFOMOD

!-----------------------------------------------------------------------
      UN = 1.d0 ; DEUX = 2.d0 ; QUATRE = DEUX*DEUX
      XLIM = 1.d-12 ; POINT5 = .5d0
      PI  = QUATRE*ATAN(UN)

      IR  = 2*NZERO ; FI0 = DBLE(IR)
      FI1 = FI0+UN  ; FN  = DBLE(NZERO)
 
      DN   = FI0/SQRT(QUATRE*FI0*FI0-UN) 
      DN1  = FI1/SQRT(QUATRE*FI1*FI1-UN) 
      A    = DN1*FI0
      B    = DN *FI1
      IRP  = IR + 1
      IRM  = IR - 1 

      SUMW16 = 0.d0 ; SUMW8 = 0.d0

      DO  I=1,NZERO 

          DI = DBLE(I) ; DOT =  DI-UN
          FI  = -PI*POINT5*(DOT+POINT5)/FN + PI*POINT5
          FI  =  SIN(FI) 

  100     CALL ordleg16 (G,FI,IR)
              CALL ordleg16 (GM,FI,IRM)
              CALL ordleg16 (GP,FI,IRP)
              GT    = (A*GP-B*GM)/(FI*FI-UN) 
              FTEMP = FI - G/GT 
              GTEMP = FI - FTEMP
              FI    = FTEMP
          IF (ABS(GTEMP) > XLIM) GOTO 100
 
          C       = DEUX*(UN-FI*FI) 
          CALL ordleg16 (D,FI,IRM) 
          D       = D*D*FI0*FI0

          C       = C*(FI0-POINT5)/D
          D       = SQRT(UN-FI*FI)

          F(I)    = FI
          WT(I)   = C
          RAD(I)  = ACOS(FI) 
          SIA(I)  = D
          WOCS(I) = C/(D*D)

          SUMW16 = SUMW16 + C ; SUMW8 = SUMW8 + WT(I)

      END DO

      CALL GET_INFOMOD( INFO )

      IF (INFO)                                    THEN
          ERRW8 = (SUMW16 - UN)
          PRINT *,' GAUSSG16: Calculation for NZERO =',NZERO
          PRINT *,'           Summ(16) Weights - 1. = ',ERRW8
          ERRW8 = (SUMW8 - 1.0)
          PRINT *,'           Summ(8)  Weights - 1. = ',ERRW8
      END IF

      RETURN
!-----------------------------------------------------------------------

      END 
      SUBROUTINE ordleg16 (SX,COA,IR)
#     if !defined (__HOS_AIX__) && !defined (__INTEL_COMPILER_UPDATE) && !defined (__GFORTRAN__)
      USE ddmodule
#     endif
      IMPLICIT NONE
 
      INTEGER        IR
#     if defined (__HOS_AIX__) || defined (__INTEL_COMPILER_UPDATE) || defined (__GFORTRAN__)
      REAL*16        SX,COA
#     else
      TYPE (DD_REAL) SX,COA
#     endif
!**    THIS ROUTINE IS A SUBSET OF BELOUSOVS ALGORITHM 
!**    USED TO CALCULATE ORDINARY LEGENDRE POLYNOMIALS
!**    AND IS CALLED BY GAUSSG16
 
!**     SX = LEGENDRE POLYNOMIAL EVALUATED AT COA
!**    COA = COSINE OF COLATITUDE
!**     IR = WAVE NUMBER 

#     if defined (__HOS_AIX__) || defined (__INTEL_COMPILER_UPDATE) || defined (__GFORTRAN__)
      REAL*16        ZERO,UN,DEUX,SQR2, PI,THETA,C1, &
                     FN,FN2,FN2SQ, ANG,S1,C4,A,B,FK, &
                     DELTA,SIA, DN
#     else
      TYPE (DD_REAL) ZERO,UN,DEUX,SQR2, PI,THETA,C1, &
                     FN,FN2,FN2SQ, ANG,S1,C4,A,B,FK, &
                     DELTA,SIA, DN
#     endif
      INTEGER        IRPP,IRPPM,N,N1,K,KK

!-----------------------------------------------------------------------
      ZERO = 0.d0 ; UN = 1.d0 ; DEUX = 2.d0 ;  SQR2  = SQRT(DEUX) 

      PI = DEUX*DEUX*ATAN(UN)
      IRPP  = IR + 1 
      IRPPM = IRPP - 1
      DELTA = ACOS(COA) 
      SIA   = SIN(DELTA) 
 
      THETA = DELTA 
      C1    = SQR2 
 
      DO  N=1,IRPPM 
          DN = DBLE(N) ; FN2 = DEUX*DN ; FN2SQ = FN2*FN2 
          C1    = C1* SQRT(UN-UN/FN2SQ)
      END DO
 
      N   = IRPPM 
      FN  = DBLE(N)
      ANG = FN*THETA
      S1  = ZERO
      C4  = UN
      A   =-UN
      B   = ZERO
      N1  = N+1

      DO  KK=1,N1,2 
          K   = KK-1
          IF (K == N) C4 = C4/DEUX
          S1  = S1+C4* COS(ANG)
          A   = A+DEUX 
          B   = B+UN
          FK  = DBLE(K)
          ANG = THETA*(FN-FK-DEUX) 
          C4  = (A*(FN-B+UN)/(B*(FN2-A)))*C4
      END DO
 
      SX = S1*C1
 
      RETURN
!-----------------------------------------------------------------------

      END 
      SUBROUTINE trigl (ILATH,SR,WR,CR,RADR,WOSQ)

      IMPLICIT none

      INTEGER   ILATH
      REAL(8)   SR(ILATH*2),  WR(ILATH*2),CR(ILATH*2), &
              RADR(ILATH*2),WOSQ(ILATH*2)

!**    JAN 19/78 - J.D.HENDERSON

!**    THE ARGUMENT LIST IS THE SAME AS FOR GAUSSG, EXCEPT THAT RADR
!**    READS IN COLATITUDE AND RETURNS THE LATITUDE. GAUSSG FILLS ONLY
!**    THE N HEM ORDERED N TO S. THIS ROUTINE MAKES THE ARRAYS GLOBAL
!**    AND ORDERED FROM S TO N. 

!**    PARAMETERS...
!**         SR   = SIN(LAT),  
!**         CR   = COS(LAT), 
!**         RADR = LATITUDE IN RADIANS.
!**         WR   = GAUSSIAN WEIGHTS, 
!**         WOSQ = WR/(SR**2).

      INTEGER ILAT, J,K
      REAL(8) PIH

!-------------------------------------------------------------------- 
!**    CR,WR,WOSQ ARE SYMMETRIC ABOUT THE EQUATOR.
!**    SR AND RAD ARE ANTISYMMETRIC.

      PIH  = 2.D0*ATAN(1.D0)

      ILAT = ILATH*2
      DO  J=1,ILATH
          K       = ILAT+1-J
          CR(K)   = CR(J)
          WR(K)   = WR(J)
          WOSQ(K) = WOSQ(J)
          SR(K)   = SR(J)
          SR(J)   =-SR(J)
          RADR(K) = PIH-RADR(J)
          RADR(J) =-RADR(K)
      END DO

      RETURN
!-------------------------------------------------------------------- 

      END 
      SUBROUTINE trigl2 (ILATH,SR,WR,CR,RADR,WOSQ,KIND)

      IMPLICIT none

      INTEGER   ILATH,KIND
      REAL(8)   SR(ILATH*2),  WR(ILATH*2),CR(ILATH*2), &
              RADR(ILATH*2),WOSQ(ILATH*2)

!**    DEC 24/92 - B.Dugas, RPN (BASED ON TRIGL BY J.D.HENDERSON(1978))

!**    THE ARGUMENT LIST IS THE SAME AS FOR TRIGL, EXCEPT FOR KIND
!**    WHICH DENOTES THE HEMISPHERIC CHARACTERISTIC OF THE DESTINATION
!**    GRID (0=GLOBAL,1=NH,2=SH). GAUSSG FILLS ONLY THE N HEM ORDERED 
!**    N TO S. THIS ROUTINE CORRECT THE ARRAYS (I.E. GLOBAL ARRAYS
!**    ARE ORDERED FROM S TO N, ETC...). 

!**    PARAMETERS...
!**         SR   = SIN(LAT),  
!**         CR   = COS(LAT), 
!**         RADR = LATITUDE IN RADIANS.
!**         WR   = GAUSSIAN WEIGHTS, 
!**         WOSQ = WR/(SR**2),
!**         KIND = HEMISPHERIC SIGNATURE.

      REAL(8) PIH,VALEUR
      INTEGER J,K,ILAT,ILATP1

!-------------------------------------------------------------------- 
      PIH  = 2.D0*ATAN(1.D0)

      IF (KIND == 0)                                           THEN

!**        GLOBAL DATA.
!**        CR,WR,WOSQ ARE SYMMETRIC ABOUT THE EQUATOR.
!**        SR AND RAD ARE ANTISYMMETRIC.


          ILAT = ILATH*2
          DO  J=1,ILATH
              K       = ILAT+1-J
              CR(K)   = CR(J)
              WR(K)   = WR(J)
              WOSQ(K) = WOSQ(J)
              SR(K)   = SR(J)
              SR(J)   =-SR(J)
              RADR(K) = PIH-RADR(J)
              RADR(J) =-RADR(K)
          END DO

      ELSE 

!**        HEMISPHERIC DATA. CHANGE THE ORDER OF THE FIELDS 
!**        AS REQUIRED. NOTE THAT THE GAUSSIAN WEIGHTS HAVE
!**        TO BE DOUBLED SINCE THERE IS ONLY ONE HEMISPHERE.
  
          ILAT = ILATH

          IF (KIND == 1)                                       THEN
  
!**            RE-ORDER FIELDS FROM EQUATOR TO NORTH POLE. 

              ILATP1 = ILAT+1
              DO  J=1, ILATP1/2

                  VALEUR         = SR(ILATP1-J)
                  SR(ILATP1-J)   = SR(       j)
                  SR(       j)   = VALEUR

                  VALEUR         = WR(ILATP1-J)*2.0
                  WR(ILATP1-J)   = WR(       j)*2.0
                  WR(       j)   = VALEUR

                  VALEUR         = CR(ILATP1-J)
                  CR(ILATP1-J)   = CR(       j)
                  CR(       j)   = VALEUR

                  VALEUR         = PIH-RADR(ILATP1-J)
                  RADR(ILATP1-J) = PIH-RADR(       j)
                  RADR(       j) = VALEUR

                  VALEUR         = WOSQ(ILATP1-J)*2.0
                  WOSQ(ILATP1-J) = WOSQ(       j)*2.0
                  WOSQ(       j) = VALEUR
              
              END DO

          ELSE IF (KIND == 2)                                  THEN
  
!**            RECONSIDER ANTI-SYMMETRIC FIELDS AS GAUSSG
!**            PRODUCES NH-CONSISTENT VALUES.  THE FIELDS
!**            ARE ALREADY WELL ORDERED (S TO N).
  
              DO  J=1,ILAT 
                  SR(J)   = -SR(J)
                  RADR(J) =  RADR(J)-PIH
                  WR(J)   =  WR(J)+WR(J)
                  WOSQ(J) =  WOSQ(J)+WOSQ(J)
              END DO
  
          END IF

      END IF

      RETURN
!-------------------------------------------------------------------- 

      END 
      SUBROUTINE gaussg (NZERO,F,WT,SIA,RAD,WOCS)
 
      IMPLICIT NONE
 
      INTEGER  NZERO
      real(8)  F(NZERO),WT(NZERO),SIA(NZERO),RAD(NZERO),WOCS(NZERO)

!**    THIS ROUTINE CALCULATES THE ROOTS (F) OF THE ORDINARY LEGENDRE
!**    POLYNOMIALS OF ORDER 2*NZERO.  THE FIRST STEP IS TO MAKE AN 
!**    INITIAL GUESS FOR EACH ROOT AND THEN TO USE THE ORDINARY
!**    LEGENDRE ALGORITHM (ordleg8) AND NEWTONS METHOD TO REFINE
!**    THE SOLUTION UNTIL THE CRITERION XLIM IS SATISFIED. THIS
!**    IS DONE ON THE FIRST HALF OF THE RANGE DUE TO THE SYMMETRY
!**    THE ZEROS WITH RESPECT TO F=0. THIS IS A REAL(8) VERSION.

!**       F = COSINE OF COLATITUDE 
!**      WT = CORRESPONDING GAUSSIAN WEIGHT
!**     SIA = SINE OF COLATITUDE 
!**     RAD = COLATITUDE IN RADIANS
!**    WOCS = GAUSSIAN WEIGHT / COS(COLAT)**2

      INTEGER   IR,IRP,IRM,I
      REAL(8)   FI,FI0,FI1,FN, G,GM,GP,GT, FTEMP,GTEMP, &
                A,B,C,D, DN,DN1, PI,XLIM, DOT, ERRW8,SUMW8

      LOGICAL   INFO
      EXTERNAL  GET_INFOMOD

!-----------------------------------------------------------------------
      If (NZERO >= 100)                                        Then
          ! Always call the extended precision
          ! version for higher truncations
          CALL gaussg16( NZERO, F,WT,SIA,RAD,WOCS )
          RETURN
      End If

      XLIM = 1.D-13 
      PI   = 4.D0*ATAN(1.D0)
      IR   = 2*NZERO
      FI0  = DBLE(IR)
      FI1  = FI0+1.D0 
      FN   = DBLE(NZERO) 
 
      DN   = FI0/SQRT(4.D0*FI0*FI0-1.D0) 
      DN1  = FI1/SQRT(4.D0*FI1*FI1-1.D0) 
      A    = DN1*FI0
      B    = DN *FI1
      IRP  = IR + 1
      IRM  = IR - 1 

      SUMW8 = 0.d0

      DO  I=1,NZERO 

          DOT =  DBLE(I-1)
          FI  = -PI*.5D0*(DOT+.5D0)/FN + PI*.5D0
          FI  =  SIN(FI) 

  100     CALL ordleg8 (G,FI,IR)
              CALL ordleg8 (GM,FI,IRM)
              CALL ordleg8 (GP,FI,IRP)
              GT    = (A*GP-B*GM)/(FI*FI-1.D0) 
              FTEMP = FI - G/GT 
              GTEMP = FI - FTEMP
              FI    = FTEMP
          IF (ABS(GTEMP) .GT. XLIM) GOTO 100
 
          C       = 2.D0*(1.D0-FI*FI) 
          CALL ordleg8 (D,FI,IRM) 
          D       = D*D*FI0*FI0

          C       = C*(FI0-.5D0)/D
          D       = SQRT(1.D0-FI*FI)

          F(I)    = FI
          WT(I)   = C
          RAD(I)  = ACOS(FI) 
          SIA(I)  = D
          WOCS(I) = C/(D*D)

          SUMW8 = SUMW8 + WT(I)

      END DO
 
      CALL GET_INFOMOD( INFO )

      IF (INFO)                                    THEN
          ERRW8 = (SUMW8 - 1.0)
          PRINT *,' GAUSSG8: Calculation for NZERO =',NZERO
          PRINT *,'          Summ(8)  Weights - 1. = ',ERRW8
      END IF

      RETURN
!-------------------------------------------------------------------- 

      END 
      SUBROUTINE ordleg8 (SX,COA,IR)
 
      IMPLICIT NONE
 
      INTEGER  IR
      REAL(8)  SX,COA

!**    THIS ROUTINE IS A SUBSET OF BELOUSOVS ALGORITHM 
!**    USED TO CALCULATE ORDINARY LEGENDRE POLYNOMIALS
!**    AND IS CALLED BY GAUSSG.
 
!**     SX = LEGENDRE POLYNOMIAL EVALUATED AT COA
!**    COA = COSINE OF COLATITUDE
!**     IR = WAVE NUMBER 

      REAL(8)  SQR2, PI,THETA,C1, FN,FN2,FN2SQ, &
               ANG,S1,C4,A,B,FK, DELTA,SIA
      INTEGER  IRPP,IRPPM,N,N1,K,KK
 
!-----------------------------------------------------------------------
      PI    = 4.D0*ATAN(1.D0)
      SQR2  = SQRT(2.D0) 
      IRPP  = IR + 1 
      IRPPM = IRPP - 1
      DELTA = ACOS(COA) 
      SIA   = SIN(DELTA) 
 
      THETA = DELTA 
      C1    = SQR2 
 
      DO  N=1,IRPPM 
          FN2   = DBLE(2*N) 
          FN2SQ = FN2*FN2 
          C1    = C1* SQRT(1.D0-1.D0/FN2SQ)
      END DO
 
      N   = IRPPM 
      FN  = DBLE(N)
      ANG = FN*THETA
      S1  = 0.D0
      C4  = 1.D0
      A   =-1.D0
      B   = 0.D0 
      N1  = N+1

      DO  KK=1,N1,2 
          K   = KK-1
          IF (K.EQ.N) C4 = 0.5D0*C4 
          S1  = S1+C4* COS(ANG)
          A   = A+2.D0 
          B   = B+1.D0 
          FK  = DBLE(K) 
          ANG = THETA*(FN-FK-2.D0) 
          C4  = (A*(FN-B+1.D0)/(B*(FN2-A)))*C4
      END DO
 
      SX = S1*C1
 
      RETURN
!-------------------------------------------------------------------- 

      END 
