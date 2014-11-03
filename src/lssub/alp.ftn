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
C     $Log: alp.ftn,v $
C     Revision 3.5  2014/10/16 12:00:35  dugas
C     Modifications commandees par la compilation avec GFORTRAN et IFORT.
C
C     Revision 3.4  2014/09/25 18:31:49  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.3  2013/10/28 20:42:42  bernard
C     Enlever du code commente/experimental.
C
C     Revision 3.2  2013/10/08 00:57:07  bernard
C      - Substiruer COMPLEX(:) par REAL(8)(2,:) dans CALPHA
C      - Ajouter les routines ALPINI2,ALPCAL2 qui sont des versions
C        adaptees des routines du SEF pour calculer les polynomes de
C        legendre selon les formules de Belousov (1962)
C      - Ajouter ALPDY3, une version de ALPDY2 adaptee
C        a cette meme formule
C
C     Revision 3.1  1994/11/17 14:12:55  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:55:12  13:55:12  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:31:31  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.0  92/02/21  11:30:47  armnrbd
C     Initial revision
C     

      SUBROUTINE alpas2 (ALP,LALP,LM,WRKS) 

C     *****   OCT 1975  -  JOHN D. HENDERSON  ****

C     * ALP(LALP,LM) CAN CONTAIN THE LEGENDRE POLYNOMIALS CALCULATED
C     *  BY SUBROUTINE ALPMN2, OR THEIR N-S DERIVATIVES CALCULATED
C     *  BY SUBROUTINE ALPDR2.

C     * THE SYMMETRIC AND ANTISYMMETRIC VALUES IN EACH ROW OF ALP 
C     * ARE ORIGINALLY INTERLEAVED. E.G. ROW 1 = (0,1,2,3,...LALP). 
C     * THIS SUBROUTINE SEPARATES EACH ROW INTO TWO PARTS.
C     *  E.G. ROW 1 = (0,2,4...LALP/2, 1,3,5...LALP ).
C     * WRKS IS AN SCM WORK FIELD OF LALP WORDS.
 
C     * WARNING - LALP MUST BE EVEN.

      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

      DIMENSION ALP(LALP,1) 
      DIMENSION WRKS(*) 

C-----------------------------------------------------------------------
      LALPH  = LALP/2
      LALPH1 = LALPH+1
 
      DO 30 M=1,LM
 
C         * TRANSFER ONE ROW OF ALP TO WRKS.
 
          DO 20 N=1,LALP
              WRKS(N) = ALP(N,M)
   20     CONTINUE
 
C         * PUT THE SYMMETRIC VALUES IN WORDS 1 TO LALP/2 OF ALP. 
 
          NSYM = -1 
          DO 22 N=1,LALPH 
              NSYM     = NSYM+2 
              ALP(N,M) = WRKS(NSYM) 
   22     CONTINUE
 
C         * PUT THE ANTISYMMETRIC VALUES IN WORDS LALP/2+1 TO LALP. 
 
          NASM = 0
          DO 24 N=LALPH1,LALP 
              NASM     = NASM+2 
              ALP(N,M) = WRKS(NASM) 
   24     CONTINUE
 
   30 CONTINUE
 
      RETURN
      END 

      SUBROUTINE alpdl2 (DELALP,ALP,LSR,LM)
 
C     * JUL 17/79 - J.D.HENDERSON 

C     * COMPUTES LAPLACIAN OF LEGENDRE POLYNOMIALS. 
C     * LSR CONTAINS ROW LENGTH INFO. 
 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

      DIMENSION DELALP(*),ALP(*), LSR(2,*)

C-------------------------------------------------------------------- 
      DO 210 M=1,LM 
          KL = LSR(2,M) 
          KR = LSR(2,M+1)-1 
          DO 210 K=KL,KR
              FNS       =  DBLE( (M-1)+(K-KL) )
              DELALP(K) = -FNS*(FNS+1.D0)*ALP(K)
  210 CONTINUE
 
      RETURN
      END 

      SUBROUTINE alpdr2 (DALP,ALP,LALP,LM,EPSI)
 
C     *****   OCT 1975  -  JOHN D. HENDERSON  ****

C     * CALCULATES N-S DERIVATIVES OF EACH ASSOCIATED LEGENDRE POLYNOMIAL

C     * DALP(LALP,LM)   WILL CONTAIN N-S DERIVATIVE OF ALP. 
C     * ALP(LALP,LM)    CONTAINS LEGENDRE POLYNOMIALS FOR ONE LATITUDE.
C     * EPSI(LALP,LM)   CONTAINS PREVIOUSLY CALCULATED CONSTANTS. 
 
C     * WARNING - LALP MUST BE EVEN.
C     *         - LAST ELEMEMT OF EACH ROW IS SET TO ZERO.
 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

      DIMENSION DALP(LALP,1),ALP(LALP,1),EPSI(LALP,1) 

C-----------------------------------------------------------------------
      LALPM = LALP-1
 
      DO 30 M=1,LM
 
          DO 20 N=1,LALPM 
              FNS    = DBLE(M+N-2)
              ALPILM = 0.D0 
              IF (N.GT.1) ALPILM = ALP(N-1,M)
              DALP(N,M) = (FNS+1.D0)*EPSI(N,  M)*ALPILM 
     +                  -  FNS      *EPSI(N+1,M)*ALP(N+1,M)
   20     CONTINUE
 
          DALP(LALP,M) = 0.D0 

   30 CONTINUE
 
      RETURN
      END 

      SUBROUTINE alpdy2 (DALP,ALP,LSR,LM,EPSI) 
 
C     * JUL 16/79 - J.D.HENDERSON 

C     * SETS DALP TO N-S DERIVATIVE OF LEGENDRE POLYNOMIALS IN ALP. 
C     * EPSI CONTAINS PRECOMPUTED CONSTANTS.
C     * LSR CONTAINS ROW LENGTH ONFO. 
 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

      DIMENSION  DALP(*),ALP(*),EPSI(*), LSR(2,*)

C-------------------------------------------------------------------- 
      DO 230 M=1,LM 
          MS = M-1
          KL = LSR(2,M) 
          KR = LSR(2,M+1)-2 
 
          DO 220 K=KL,KR
              FNS    = DBLE(MS+K-KL)
              ALPILM = 0.D0 
              IF (K.GT.KL) ALPILM = ALP(K-1) 
              DALP(K) = (FNS+1)*EPSI(K)  *ALPILM
     +                -  FNS   *EPSI(K+1)*ALP(K+1)
  220     CONTINUE

          DALP(KR+1) = 0.D0

  230 CONTINUE
 
      RETURN
      END 

      SUBROUTINE alpn2 (ALP,NLP,ILAT,SL,EPSI)
 
C     * JULY 10/78 - J.D.HENDERSON

C     * COMPUTES LEGENDRE POLYNOMIALS.
 
C     * ALP(NLP,ILAT) = LEGENDRE POLYNOMIALS (NLP VALUES AT ILAT LATS). 
C     * SL(ILAT)      = SINE OF THE LATITUDES. 
 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

      DIMENSION ALP(NLP,ILAT),SL(ILAT),EPSI(NLP) 

C-------------------------------------------------------------------- 
      EPSI(1)      = 0.D0
      DO 110 N=2,NLP
          NS       = N-1
          EPSI(N)  = NS/SQRT(4.D0*NS*NS-1.D0)
  110 CONTINUE
 
      DO 210 J=1,ILAT 
          SINLAT   = SL(J)
          ALP(1,J) = 1.D0/SQRT(2.D0)
          ALP(2,J) = SQRT(3.D0/2.D0)*SINLAT 
          DO 210 N = 3,NLP
              ALP(N,J) = ( SINLAT*ALP(N-1,J) - EPSI(N-1)*ALP(N-2,J) )
     +                 /                  EPSI(N) 
  210 CONTINUE
 
      RETURN
      END 

      SUBROUTINE alpnd2 (DALP,ALP,NLP,ILAT,EPSI) 
 
C     * JUL 12/78 - J.D.HENDERSON 

C     * COMPUTES DERIVATIVES OF LEGENDRE POLYNOMIALS. 
 
C     * ALP(NLP,ILAT)  = LEGENDRE POLYNOMIALS (NLP VALUES AT ILAT LATS). 
C     * DALP(NLP,ILAT) = CORRESPONDING DERIVATIVES OF THE POLYNOMIALS.

      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
 
      DIMENSION DALP(NLP,ILAT),ALP(NLP,ILAT), EPSI(NLP)

C-------------------------------------------------------------------- 
C     * PRECOMPUTE CONSTANTS IN EPSI. 
 
      NLPM        = NLP-1

      EPSI(1)     = 0.D0
      DO 110 N=2,NLP
          NS      = N-1
          EPSI(N) = NS/SQRT(4.D0*NS*NS-1.D0)
  110 CONTINUE
 
      DO 210 J=1,ILAT 
          DALP(1,J)     = 0.D0
          DALP(NLP,J)   = 0.D0
          DO 210 N=2,NLPM 
              FNS       = DBLE(N-1)
              DALP(N,J) = (FNS+1.D0)*EPSI(N)  *ALP(N-1,J)
     +                  -  FNS      *EPSI(N+1)*ALP(N+1,J)
  210 CONTINUE
 
      RETURN
      END 

      SUBROUTINE alpnm2 (ALP,LALP,LM,SINLAT,EPSI)
 
C     *****   OCT 1975  -  JOHN D. HENDERSON  ****

C     * PUTS LEGENDRE POLYNOMIALS IN ALP(LALP,LM)  FOR ONE LATITUDE. 
C     * SINLAT IS THE SINE OF THE REQUIRED LATITUDE. EPSI IS A FIELD
C     * OF CONSTANTS THE SAME SIZE AS ALP.
 
C     * THE SYMMETRIC AND ANTISYMMETRIC VALUES IN ALP ARE INTERLEAVED 
C     * IN EACH ROW.   I.E.  ROW 1 = ( 0,1,2,3....LALP ). 
 
C     * WARNING - LALP MUST BE EVEN AND ITS MINIMUM VALUE IS 4.
 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

      DIMENSION ALP(LALP,1),EPSI(LALP,1)

C-----------------------------------------------------------------------
      COS2 = 1.D0 - SINLAT*SINLAT 
      PROD = 1.D0 
      A    = 1.D0
      B    = 0.D0
 
C     * LOOP 30 COVERS LONGITUDINAL WAVE NUMBERS 0 TO LM-1. 
 
      DO 30 M=1,LM
          FM = DBLE(M-1) 
          IF (M.NE.1)                                          THEN
              A    = A+2.D0
              B    = B+2.D0
              PROD = PROD*COS2*A/B
          END IF
 
C         * COMPUTE THE FIRST TWO ELEMEMTS OF THE ROW.
 
          ALP(1,M) = SQRT(0.5D0*PROD)
          ALP(2,M) = SQRT(2.D0*FM+3.D0)*SINLAT*ALP(1,M) 
 
C         * NOW COMPUTE ELEMENTS 3 TO LR IN THE ROW IN PAIRS. 
 
          DO 20 N=3,LALP,2
              ALP(N,  M) = (SINLAT*ALP(N-1,M) - EPSI(N-1,M)*ALP(N-2,M)) 
     +                   /                EPSI(N,M) 
              ALP(N+1,M) = (SINLAT*ALP(N,  M) - EPSI(N,  M)*ALP(N-1,M))
     +                   /               EPSI(N+1,M) 
   20     CONTINUE
 
   30 CONTINUE
 
      RETURN
      END 

      SUBROUTINE alpst2 (ALP,LSR,LM,SINLAT,EPSI) 
 
C     * JUL  2/79 - J.D.HENDERSON 

C     * PUTS  LEGENDRE POLYNOMIALS IN ALP FOR ONE LATITUDE.  SINLAT
C     * IS THE SINE OF THE LATITUDE. EPSI IS A FIELD OF PRECOMPUTED
C     * CONSTANTS.  LSR CONTAINS ROW LENGTH INFO FOR ALP,EPSI.
 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

      DIMENSION LSR (2,*), ALP(*),EPSI(*) 

C-------------------------------------------------------------------- 
      COS2 = 1.D0-SINLAT*SINLAT 
      PROD = 1.D0
      A    = 1.D0
      B    = 0.D0
 
C     * LOOP OVER LONGITUDINAL WAVE NUMBERS.
 
      DO 230 M=1,LM 
          FM = DBLE(M-1) 
          IF (M.NE.1) THEN
              A    = A+2.D0
              B    = B+2.D0
              PROD = PROD*COS2*A/B
          END IF
 
C         * COMPUTE THE FIRST VALUE IN THE ROW AND
C         * COMPUTE THE SECOND ELEMENT ONLY IF NEEDED.
 
          KL = LSR(2,M) 
          KR = LSR(2,M+1)-1 

          ALP(KL) = SQRT(0.5D0*PROD) 
          IF (KR.GT.KL) ALP(KL+1) = SQRT(2.D0*FM+3.D0)*SINLAT*ALP(KL)
 
C         * COMPUTE THE REST OF THE VALUES IN THE ROW IF NEEDED.
 
          KL2 = KL+2
          IF (KL2.LE.KR)                                       THEN
              DO 220 K=KL2,KR 
                  ALP(K) = ( SINLAT*ALP(K-1) - EPSI(K-1)*ALP(K-2) )
     +                   /                EPSI(K)
  220         CONTINUE
          END IF
 
  230 CONTINUE

      RETURN
      END 

      SUBROUTINE alpwei (ALP,DALP,DELALP, WL,WOSL,LSR,LM) 
 
C     * JUL 21/83 - R.LAPRISE.

C     * NORMALIZE  THE LEGENDRE POLYNOMIALS AND THEIR DERIVATIVES
C     * BY THE GAUSSIAN WEIGTHS AS REQUIRED FOR THE DIRECT TRANS-
C     * FORM FROM FOURIER TO SPECTRAL.  DELALP ALSO INCLUDES  THE 
C     * FACTOR -0.5 REQUIRED FOR THE DIVERGENCE TENDENCY TERM IN-
C     * VOLVING THE KINETIC ENERGY.
 
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

      DIMENSION  ALP(*),DALP(*),DELALP(*), LSR(2,*)
C-----------------------------------------------------------------------
      NR = LSR(2,LM+1)-1
      DO 100 I=1,NR 
           ALP(I) =    ALP(I)*WL
          DALP(I) =   DALP(I)*WOSL
        DELALP(I) = DELALP(I)*WOSL*(-0.5D0) 
  100 CONTINUE

      RETURN
      END 

      SUBROUTINE calpha (ALPHA, LSR,LM,LA, KIND)
 
C     * NOV 26/90 - B.DUGAS, RPN.

C     * Calpha RETOURNE CMPLX-ALPHA=(NS*(NS+1),0), TOUT EN
C     * CONSIDERANT LA POSSIBLILITE DE CHAMPS SYMMETRIQUES
C     * ( SI KIND>0 ) OU ANTI-SYMMETRIQUES  (SI KIND<0 ).

      IMPLICIT none

      INTEGER  LSR(2,*),LM,LA,KIND
      REAL(8)  ALPHA(2,LA) 

      INTEGER  PAS,DEB,M,N,NI,NF,L

C-----------------------------------------------------------------------
C     * CONSIDERER LES CORRECTIONS DUES AUX CAS KIND.NE.0

          DEB = 0
      IF (KIND.EQ.0)                                           THEN
          PAS = 1
      ELSE
          PAS = 2
          IF (KIND.LT.0) DEB = 1
      END IF

C     * TRAITER M=0.

      L  = 0
      NI = LSR(2,1)
      NF = LSR(2,2)-2
      DO 300 N=0,NF-NI
          L     = L+1
          ALPHA(1,L) = N*(N+1) ; ALPHA(2,L) = 0.0
  300 CONTINUE

C     * TRAITER M.NE.0

      DO 350 M=1,LM-1
          NI = LSR(2,M+1)+DEB
          NF = LSR(2,M+2)-2
          DO 350 N=0,NF-NI,PAS
              L     = L+1
              ALPHA(1,L) = (M+N)*(M+N+1) ; ALPHA(2,L) = 0.0
  350 CONTINUE

C-----------------------------------------------------------------------
      RETURN
      END 

      SUBROUTINE epscal (EPSI,LSR,LM)

C     * JUL  2/79 - J.D.HENDERSON 

C     * COMPUTES ARRAY EPSILON USED IN THE SPECTRAL MODEL.
C     * EPSILON(N,M)=SQRT((N**2-M**2)/(4*N**2-1)) 
C     * LSR CONTAINS ROW START INFO.

      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

      DIMENSION  EPSI(*), LSR(2,*)

C---------------------------------------------------------------------
      DO 220 M=1,LM 
          MS = M-1
          KL = LSR(2,M) 
          KR = LSR(2,M+1)-1 

          K1 = KL 
          IF (M.EQ.1) K1=2 

          DO 220 K=K1,KR
              NS      = MS+(K-KL)
              FNUM    = DBLE(NS*NS-MS*MS) 
              FDEN    = DBLE(4*NS*NS-1) 
              EPSI(K) = SQRT(FNUM/FDEN) 
  220 CONTINUE

      EPSI(1) = 0.D0

      RETURN
      END 

      SUBROUTINE epsil2 (EPSI,LALP,LM) 

C     *****   OCT 1975  -  JOHN D. HENDERSON  ****

C     * CALCULATES EPSILON(N,M) = SQRT((N**2-M**2)/(4*N**2 - 1))
C     * FOR N FROM 0 TO LALP-1, AND M FROM 0 TO LM-1 IN EPSI(LALP,LM).
C     * EPSILON IS A FIELD OF CONSTANTS USED IN THE SPECTRAL MODELS.

      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

      DIMENSION EPSI(LALP,1)
C-----------------------------------------------------------------------
C 
      DO 20 M=1,LM
          MS = M-1
          N1 = 1
          IF (MS.EQ.0) N1 = 2

          DO 20 N=N1,LALP 
              NS        = MS+N-1 
              FNUM      = DBLE(NS*NS - MS*MS) 
              FDEN      = DBLE(4*NS*NS - 1) 
              EPSI(N,M) = SQRT(FNUM/FDEN) 
   20 CONTINUE

      EPSI(1,1)=0.D0

      RETURN
      END 

      SUBROUTINE GALP(AN,NLP,A,W,ILAT,ALP)

C     * JULY 10/78 - J.D.HENDERSON
C     * COMPUTES LEGENDRE POLYNOMIAL SERIES FROM DATA VALUES.

C     * AN(NLP)       = AMPLITUDES FOR LEGENDRE POLYNOMIAL SERIES.
C     * A(ILAT)       = DATA VECTOR DEFINED ON GAUSSIAN LATITUDES.
C     * W(ILAT)       = GAUSSIAN WEIGHTS FOR EACH GAUSSIAN LATITUDE.
C     * ALP(NLP,ILAT) = LEGENDRE POLYNOMIALS (NLP VALUES AT ILAT LATS).

      IMPLICIT none

      INTEGER  NLP,ILAT,N,J
      REAL*8   ALP(NLP,ILAT),SUM
      REAL*8   AN(NLP),A(ILAT),W(ILAT)
C--------------------------------------------------------------------

      DO 300 N=1,NLP
          SUM = 0.
          DO 200 J=1,ILAT
              SUM = SUM+W(J)*A(J)*ALP(N,J)
  200     CONTINUE
          AN(N) = SUM
  300 CONTINUE

      RETURN
      END

      SUBROUTINE LPAG(A,ILAT,AN,NLP,ALP)

C     * JULY 10/78 - J.D.HENDERSON
C     * COMPUTES DATA VALUES FROM LEGENDRE POLYNOMIAL SERIES.

C     * AN(NLP)       = AMPLITUDES FOR LEGENDRE POLYNOMIAL SERIES.
C     * A(ILAT)       = DATA VECTOR DEFINED ON GAUSSIAN LATITUDES.
C     * ALP(NLP,ILAT) = LEGENDRE POLYNOMIALS (NLP VALUES AT ILAT LATS).

      IMPLICIT none

      INTEGER  NLP,ILAT,N,J
      REAL*8   ALP(NLP,ILAT),SUM
      REAL*8   AN(NLP),A(ILAT)
C--------------------------------------------------------------------

      DO 300 J=1,ILAT
          SUM = 0.
          DO 200 N=1,NLP
              SUM = SUM+AN(N)*ALP(N,J)
  200     CONTINUE
          A(J) = SUM
  300 CONTINUE

      RETURN
      END
      SUBROUTINE ALPINI2 (LSR,LM)

      IMPLICIT none

      INTEGER  LM,LSR(2,LM+1)

!**    CALCULATES EPSILON AND OTHER CONSTANTS NEEDED IN LEGENDRE
!**    POLYNOMIAL RECURSIVE CALCULATIONS FOR N FROM 0 TO LALP(M)-1,
!**    AND M FROM 0 TO LM-1.

!**    Les resultats sont places dans le bloc /ALPI/.

!**    ****   OCT 1975  -  JOHN D. HENDERSON  ****

!**    B.DUGAS    janvier 1994
!**               Calcul des constantes c,d,e utilisees par la
!**               formule de Belousov pour calculer les ALP.
!**               Voir les details dans MWR[118], PP. 2248-2251.
!**    B.Dugas    octobre 2013
!**               Version adaptee a r.diag

      REAL*8    FNUM,FDEN,c,d,e
      INTEGER   M,N,N1,LR,MS,NS,IER, NALP0

#     include  "calpi.cdk"

      c(m,n) = sqrt( (dble( 2*n+1 )*(n+m-1)*(n+m-3)) /
     +               (dble( 2*n-3 )*(n+m  )*(n+m-2)) )
      d(m,n) = sqrt( (dble( 2*n+1 )*(n+m-1)*(n-m+1)) /
     +               (dble( 2*n-1 )*(n+m  )*(n+m-2)) )
      e(m,n) = sqrt(  dble( (2*n+1)*(n-m) )/
     +                dble( (2*n-1)*(n+m) ))

!-----------------------------------------------------------
!**   Allouer l'espace pour les AC,AD,AE, ALP, EPSI et DEBM.

      NALP0 = LSR(2,2)-1 ! Nombre de composantes en N pour M=0

      if (associated( AC )) then
          if (LSR(1,LM+1) /= LASTLA) THEN
              deallocate( AC,AD,AE,ALP,EPSI )
              if (associated( DALP )) deallocate( DALP )
              LASTLA = 0
          else
              return
          end if
      end if

      LASTLA = LSR(1,LM+1) ! Sauver pour prochain appel

      allocate( AC(NALP0,LM),AD(NALP0,LM),AE(NALP0,LM), \
                ALP(NALP0,LM),EPSI(NALP0,LM) )

      DO  M=1,LM
          MS = M-1
          N1 = 1
          IF (MS.EQ.0) N1 = 2
          LR = LSR(2,M+1)-LSR(2,M)
          DO  N=N1,LR 
              NS   = MS+N-1 
              FNUM =   NS**2 - MS**2
              FDEN = 4*NS**2 - 1 
              EPSI(N,M) = SQRT( FNUM/FDEN )
          END DO
      END DO

      EPSI(1,1) = 0.
      
!**    Remplir AC,AD et AE avec les C(M,N),D(M,N) et E(M,N)
!**    sans les calculer pour M=0,1.

      DO  M=3,LM
          DO  N=1,NALP0
              AC(N,M) = c(M-1,M-1+N-1)
              AD(N,M) = d(M-1,M-1+N-1)
              AE(N,M) = e(M-1,M-1+N-1)
          END DO
      END DO
!-----------------------------------------------------------

      RETURN
      END 
      SUBROUTINE ALPCAL2 (LSR,LM,SINLAT) 

      IMPLICIT  none

      INTEGER  LM,LSR(2,LM+1)
      REAL(8)  SINLAT
     
!**    AUTHOR: B.DUGAS, SEPT 1993
!**    (THIS IS A RECODED VERSION BASED ON 1) A FIRST VERSION 
!**     BY J.HENDERSON AND ON 2) COMMENTS REGARDING BELOUSOV'S
!**     FORMULA AS FOUND IN MWR[118], PP. 2248-2251)

!**    B.Dugas    octobre 2013
!**               Version adaptee a r.diag (NON VECTORISEE)

      INTEGER   M,N,LR
      REAL(8),  PARAMETER :: POINT5=.5,POINT75=.75

#     include  "calpi.cdk"

!---------------------------------------------------------------
      if (.not.associated( AC ) .or. LSR(1,LM+1) /= LASTLA) \
          call ALPINI2( LSR,LM )

      ! DETERMINE THE TYPE OT TRUNCATION.

      DO  M=1,LM

          LR = LSR(2,M+1)-LSR(2,M)

          IF (M == 1 .OR. M == 2) THEN

              ! COMPUTE THE FIRST ELEMENT OF THE FIRST TWO COLUMNS.

              IF (M.EQ.1) THEN
                  ALP(1,M) = SQRT( POINT5 )
              ELSE IF(M.EQ.2) THEN
                  ALP(1,M) = SQRT( POINT75*(1.-SINLAT**2) )
              ENDIF

              ! COMPUTE THE SECOND ELEMENT OF THE FIRST TWO COLUMNS.

              ALP(2,M) = (SINLAT*ALP(1,M))/EPSI(2,M)

              ! COMPUTE THE REMAINING ELEMENTS OF THE FIRST TWO COLUMNS.

              DO  N=3,LR
                  ALP(N,M)=(SINLAT*ALP(N-1,M)-EPSI(N-1,M)*ALP(N-2,M)) \
                          /              EPSI(N,M)
              END DO

          ELSE

              ! CALCULATE THE REST USING BELOUSOV'S FORMULA.

              ALP(1,M) = AC(1,M)*ALP(1,M-2)-AD(1,M)*ALP(2,M-2)*SINLAT

              DO  N=2,LR-1
                  ALP(N,M) = AC(N,M)*ALP(N  ,M-2) \
                            -AD(N,M)*ALP(N+1,M-2)*SINLAT \
                            +AE(N,M)*ALP(N-1,M  )*SINLAT
              END DO

              ! ALP(LR,M) = AC(LR,M)*ALP(LR,M-2)+AE(LR,M)*ALP(LR-1,M)
              
              ! CALCULATE THE LAST ROW USING THE OLD RECURSIVE FORMULA
                ALP(LR,M) = \
                ( SINLAT*ALP(LR-1,M) - EPSI(LR-1,M)*ALP(LR-2,M) ) \
                /                 EPSI(LR,M)
                 

          END IF

      END DO
!-----------------------------------------------------------

      RETURN
      END 
      SUBROUTINE alpdy3 (LSR,LM) 
 
      IMPLICIT  NONE

      INTEGER  LM,LSR(2,LM+1)

!     JULY 16/79 - J.D.HENDERSON : Original version (ALPDY2)
!     Octobre 2013 - B.Dugas : Version Aadaptee a Belousov pour r.diag

!     SETS DALP TO N-S DERIVATIVE OF LEGENDRE POLYNOMIALS IN ALP. 
!     EPSI CONTAINS PRECOMPUTED CONSTANTS.
!     LSR CONTAINS ROW LENGTH ONFO. 
 
      INTEGER   M,N,NALP0,LR
      REAL(8)   ALPILM,FNS

#     include  "calpi.cdk"
      
!-------------------------------------------------------------------- 
      NALP0 = LSR(2,2)-1 ! Nombre de composantes en N pour M=0

      if (.NOT.associated( DALP )) allocate( DALP(NALP0,LM) )

      DO  M=1,LM 

          LR = LSR(2,M+1)-LSR(2,M)

          DO  N=1,LR-1
              FNS    = DBLE(M+N-2)
              ALPILM = 0.D0 
              IF (N > 1) ALPILM = ALP(N-1,M) 
              DALP(N,M) = (FNS+1)*EPSI(N  ,M)*ALPILM \
                        -  FNS   *EPSI(N+1,M)*ALP(N+1,M)
          END DO

          DALP(LR,M) = 0.D0

      END DO
!-----------------------------------------------------------

      RETURN
      END 
