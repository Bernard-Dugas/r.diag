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
C     $Log: fftini.ftn,v $
C     Revision 3.8  2014/10/16 12:00:37  dugas
C     Modifications commandees par la compilation avec GFORTRAN et IFORT.
C
C     Revision 3.7  2014/09/25 18:42:02  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.6  2008/04/28 21:40:20  dugas
C     Corriger l'usage des macros pour r.gppf (passe 2).
C
C     Revision 3.5  2007/08/29 20:14:47  dugas
C     Versions F90 des routines d'initialisations.
C
C     Revision 3.4  1998/07/23 19:46:28  armnrbd
C     Corriger le calcul de TRIG dans FFTINI2.
C
C     Revision 3.3  1995/11/01  20:04:03  armnrbd
C     De-allouer les appels multi-taches.
C
C     Revision 3.2  95/02/20  20:45:23  armnrbd
C     Ajouter 45 nombres premiers a la routine FFTini2.
C     
C     Revision 3.1  94/11/17  14:13:17  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C     
C     Revision 3.0  94/11/17  13:55:27  13:55:27  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:31:41  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.5  93/08/23  11:56:04  armnrbd
C     Modifications cosmetiques.
C     
C     Revision 1.4  93/02/16  13:38:38  armnrbd
C     Modifier le comportement pour de petites valeurs de n.
C     Ajouter des nombres premiers. Corriger les commentaires.
C     
C     Revision 1.3  92/12/21  14:46:52  armnrbd
C     BugFix fonctionnel du niveau 1.2.
C     
C     Revision 1.2  92/12/21  10:41:09  armnrbd
C     Ajouter la routine FFTINI2.
C     Corriger la routine FTSET.
C     
C     Revision 1.1  92/12/17  22:07:43  armnrbd
C     Tout calculer en REAL*8.
C     
C     Revision 1.0  92/02/21  11:32:39  armnrbd
C     Initial revision
C     
#     if !defined (nombre_de_facteurs)
#         define   nombre_de_facteurs 19
#     endif
*     SUBROUTINE FFTini2 COMPUTES FACTORS OF N FOR THE STOCKHAM 
*     TRANSFORM. IT IS BASED ON 'SET77' FOUND IN RMNLIB/CRAY.
*     A MAXIMUM TOTAL OF MAXFAC FACTORS ARE CALCULATED. CONSIDERS
*     THE FIRST 30 PRIME NUMBERS AS WELL AS 16,12,10,9,8,6 AND 4.
*     IT ALSO CALCULATES THE REQUIERED TRIG FACTORS (SIN, COS AND
*     TAN [ I * (2*PI/N) ], I=0,N-1).
          
      SUBROUTINE FFTini2 (OTRIG,OFAC,M,IER) 

!!!    AUTHOR: B.Dugas - November 6 1991.

      IMPLICIT none

!!!    PARAMETRES:
!!!    OTRIG (OUT) - TRIGONOMETRIC FACTORS (REAL DIMENSION: 3 * N)
!!!    OFAC  (OUT) - N DECOMPOSITION FACTORS (REAL DIMENSION: MAXFAC+1)
!!!    N     (IN)  - TRANSFORM LENGTH. CAN BE NEGATIVE, SEE NOTE BELOW.
!!!    IER   (OUT) - ERROR FLAG, IER.NE.0 INDICATES A PROBLEM

      INTEGER    OFAC(*),M,IER
      REAL*8     OTRIG(*)

!     The output arguments OTRIG and OFAX are ignored if FFTini2 is
!     invoked with the negative of the FFT length.  The underlying
!     subroutines do not use them.

!     STKFFT8 common variables

      INTEGER    MAXFAC
      PARAMETER( MAXFAC = nombre_de_facteurs )

      CHARACTER  STRING*8
      INTEGER    N, IFAC
      REAL*8,    POINTER :: TRIG(:,:)

      COMMON   / FTSCAL  / STRING
      COMMON   / STKFFT8 / TRIG, IFAC(MAXFAC+1), N

      LOGICAL              INFO
      COMMON   /ZZVERBO  / INFO

!     Local variables

      INTEGER    JFAC(MAXFAC),LFAC(72),N0,
     +           I,K,L,IFAX,NF,NU,HOLD,ERR
      REAL*8     TPISN

      EXTERNAL   NOPRIM

      DATA       LFAC / 373,367,359,353,349,347,337,331,317, 
     +                  313,311,307,292,282,281,277,271,269, 
     +                  263,257,251,241,239,233,229,227,223, 
     +                  211,199,197,193,191,181,179,173,167, 
     +                  163,157,151,149,139,137,131,127,113, 
     +                  109,107,103,101, 97, 89, 83, 79, 73, 
     +                   71, 67, 61, 59, 53, 47, 43, 41, 37, 
     +                   31, 29, 23, 19, 17, 13, 11,  7,  1 /

!------------------------------------------------------------------ 
      N0 = ABS( M )

!!!    INITIALISE IER, ALLOCATE TRIG.

      IER = 1

      IF (string.EQ.'InitialS')                                THEN
         IF (N0.NE.N)                                          THEN
            N = N0
            DEALLOCATE( TRIG )
            ALLOCATE( TRIG(0:N-1,3) )
         ELSE
            IER = 0
            IF (M.GT.0)                                        THEN
               DO I=0,3*N-1
                  OTRIG(I+1) = TRIG(I,1)
               END DO
               DO I=1,MAXFAC+1
                  OFAC(I) = IFAC(I)
               END DO
            END IF
            RETURN
         END IF         
      ELSE
         N = N0
         ALLOCATE( TRIG(0:N-1,3) )
         string = 'InitialS'
      END IF

      IF (N.LE.12 .OR. N.EQ.16)                                THEN
         JFAC(1) = N
         NU      = 1
         K       = 1
         GOTO 030
      END IF

!!!    FIRST, FIND PRIME FACTORS FACTORS OF N. 

      NU   = N
      K    = 0 
      L    = 1 
      IFAX = LFAC(L)

  010 CONTINUE

          IF (MOD(NU,IFAX).EQ.0)                               THEN 

              K = K+1 

              IF (K.GT.MAXFAC)                                 THEN

!!!                TOO MANY FACTORS.

                  IF (INFO) WRITE(*,1000) MAXFAC
                  IFAC(1) = -1

                  RETURN

              END IF

              JFAC(K) = IFAX
              NU      = NU/IFAX

              IF (NU.EQ.1) GOTO 030 
                           GOTO 010

          END IF

  020     IFAX = LFAC(L+1)

          L    = L+1 

      IF (IFAX.GT.1) GOTO 010 

!!!    NOW CONSIDER NON PRIME FACTORS.

      CALL NOPRIM( NU,JFAC,K,MAXFAC )

      IF (NU.NE.1 .OR. K.GE.MAXFAC)                            THEN

!!!        ILLEGAL FACTORS.

          IF (INFO .AND. MAXFAC.LT.K) WRITE(*,1000) MAXFAC
          IF (INFO .AND. NU    .NE.1) WRITE(*,1001) NU
          IFAC(1) = -1

          RETURN

       END IF
!-----------------------------------------------------------------

!!!    NOW SORT FACTORS FROM SMALLEST ONWARD
!!!    AND STORE NF AS IFAC(1).

  030 CONTINUE

      IER = 0
      NF  = K

      DO I=1,NF-1
         DO L=I+1,NF
            IF (JFAC(I).GT.JFAC(L))                            THEN
                HOLD    = JFAC(L)
                JFAC(L) = JFAC(I)
                JFAC(I) = HOLD
            END IF
         END DO
      END DO

      IFAC(1) = NF
      DO I=1,NF
         IFAC(I+1) = JFAC(I)
      END DO

      IF (M.GT.0)                                              THEN
         DO I=1,MAXFAC+1
            OFAC(I) = IFAC(I)
         END DO
      END IF

      IF (INFO) WRITE(*,1002) (IFAC(I),I=1,NF+1)

!!!    INITIALISE TRIG FACTORS.

      TPISN = 4.0*ASIN(1.0D0)/N

      DO I=0,N-1
         TRIG(I,1) = COS( I*TPISN )
         TRIG(I,2) = SIN( I*TPISN )
         IF (4*I.NE.N .AND. 4*I.NE.3*N) 
     +   TRIG(I,3) = TRIG(I,2)/TRIG(I,1)
      END DO

      IF (M.GT.0)                                              THEN
         DO I=0,3*N-1
            OTRIG(I+1) = TRIG(I,1)
         END DO
      END IF

!!!    DEFINIR LA VARIABLE STRING POUR SFFT2.

      STRING = 'InitialS'

      RETURN

!-------------------------------------------------------------------
 1000 FORMAT(' Number of factors is larger than MAXFAC=',I4)
 1001 FORMAT(' NU=',I5,' contains an illegal factor.')
 1002 FORMAT(' STOCKHAM  NFAX=',I4,', factors ...',10I4)

      END

      SUBROUTINE FFTini (N,IFAC,NP,TRIG)

***    MODIFIEE LE 31 MAI 1991 PAR B.Dugas, RPN: 
***    ...  DIMENSIONNER TRIG (2*N+4*N) PLUTOT QUE (2*N*N+4*N).

***    INITIALIZE TAU AND OMEGA USED BY THE SFFT STOCKHAM.
***    THE FACTOR DECOMPOSITION IS PROVIDED BY THE ROUTINE FTSET.

      IMPLICIT  none

      INTEGER   UN,TAU,OMEG,IFAC(20),N,NP,NF,
     +          L1,L2,IP,IDO,IPW,IDW, K0,IER
      REAL*8    TRIG(6*N),PI,TPI,ARGT,ARGO
      CHARACTER STRING*8

      COMMON   /FTSCAL/ STRING

      EXTERNAL  FTset,PASini

      DATA      UN  / 1 /

*--------------------------------------------------------------------
      CALL FTset (IFAC,N,NP,IER)
*     IF (IER.NE.0) CALL                           XIT (' FFTini ',-1)

***    DEFINE POINTERS TO TAU AND OMEGA.

      TAU  = 0
      OMEG = 2*N

***    START CALC.

      NF  = IFAC(UN)
      PI  = 4.0*ATAN(1.D0)
      TPI = 2.0*PI

      L1  = UN
      IPW = UN
      IDW = UN

      DO 100 K0=UN,NF

         IP   = IFAC(K0+UN)
         ARGT = TPI/FLOAT(IP)

         L2   = IP*L1
         IDO  = N/L2
         ARGO = TPI/FLOAT(IP*IDO)

         CALL PASini (IDO,IP, TRIG(TAU+IPW),TRIG(OMEG+IDW), ARGO,ARGT)

         L1   = L2
         IPW  = IPW+IP*2
         IDW  = IDW+IDO*IP*2

  100 CONTINUE

***    YES, FFTini HAS BEEN CALLED SUCCESSFULLY.

      STRING = 'InitialS'

      RETURN
      END

      subroutine setfft8_b( n )

      implicit none

      integer  n

      external FFTini2

      real*8   prs1(1)
      integer  prs2(1),ier

      if (mod( n,2 ).ne.0) then
          print *,'Transform lenght has to be even, n=',n
          call                                     xit(' SetFFT8',-1 )
      end if

      CALL FFTini2( prs1,prs2, - abs( N ) / 2, ier ) 

      return

      end
      SUBROUTINE PASini (IDO,IP, TAU,OMEG, ARGO,ARGT)

***    MODIFIEE LE 31 MAI 1991 PAR B.Dugas, RPN: 
***    ...  DIMENSIONNER TAU(2,IP) PLUTOT QUE TAU(2,IP,IP).

***    SUBROUTINE NEEDED BY SFFTini

      IMPLICIT none

      INTEGER  IDO,IP,    I,J,L, IJ
      REAL*8   TAU(2,IP), OMEG(2,IDO,IP),
     +         ARGO,ARGT

*------------------------------------------------------------------------


      DO 6 L=1,IP

            TAU(1,L) = COS( (L-1)*ARGT )
            TAU(2,L) = SIN( (L-1)*ARGT )

    6 CONTINUE

      DO 7 J=1,IP
         DO 7 I=1,IDO
            IJ =(I-1)*(J-1)

            OMEG(1,I,J) = COS( IJ*ARGO )
            OMEG(2,I,J) = SIN( IJ*ARGO )

    7 CONTINUE

      RETURN
      END

*     SUBROUTINE FTset COMPUTES FACTORS OF N FOR THE STOCKHAM 
*     TRANSFORM. IT IS BASED ON 'SET77' FOUND IN RMNLIB/CRAY.
*     A MAXIMUM TOTAL OF MAXFAC FACTORS ARE CALCULATED. CONSIDERS
*     THE FIRST nombre_de_facteurs PRIME NUMBERS AS WELL AS 16,12,10,9,8,6 AND 4.
*     TRIGONOMETRIC FACTORS ARE CALCULATED BY FFTINI, WHICH CALLS
*     FTset.
          
      SUBROUTINE FTset (IFAC,N,NP,IER) 
  
***    AUTEUR: B.Dugas - 6 novembre 1991.
  
      IMPLICIT none

***    PARAMETRES:
***    IFAX (OUT) - VECTEUR DE FACTEURS DE DECOMPOSITION DE N.
***    N    (IN)  - NOMBRE DE POINTS DANS LA TRANSFORMEE.
***    NP   (OUT) - NOMBRE DE PROCESS PARRALLE PERMIS.
***    IER  (OUT) - DRAPEAU D'ERREUR. SI TOUT VA BIEN IL CONTIENT
***                 ZERO ET S'IL Y A UN PROBLEME, IL CONTIENT -1.

      INTEGER    MAXFAC
      PARAMETER( MAXFAC = nombre_de_facteurs )

      INTEGER  IFAC(MAXFAC+1),JFAC(MAXFAC),LFAC(16),
     +         I,K,L,IFAX,N,NF,NU,HOLD,    NP,IER

      INTEGER  mpserv, Bidon
      EXTERNAL mpserv
      EXTERNAL NOPRIM

      DATA     LFAC / 61,59,53,47,43,41,37,31,
     +                29,23,19,17,13,11,7,1 /

*------------------------------------------------------------------ 
***    CHECK THE NUMBER OF PARALLEL THREADS AVAILABLE.

                   NP    = mpserv( 'THREADS', NP   )
CC$   IF (NP.GT.1) Bidon = mpserv( 'BLOCK',Bidon )
                   NP    = MAX( NP,1 )

***    INITIALISE IER.

      IER = -1

      IF (N.LE.12 .OR. N.EQ.16)                                THEN

***       CORRECTION FOR SMALL N.

         IF (N.EQ.16)                                          THEN

             JFAC(1) = 4
             JFAC(2) = 4

         ELSE IF (N.EQ.12)                                     THEN

             JFAC(1) = 3
             JFAC(2) = 4

         ELSE IF (N.EQ.10)                                     THEN

             JFAC(1) = 2
             JFAC(2) = 5

         ELSE IF (N.EQ.9)                                      THEN

             JFAC(1) = 3
             JFAC(2) = 3

         ELSE IF (N.EQ.8)                                      THEN

             JFAC(1) = 2
             JFAC(2) = 4

         ELSE IF (N.EQ.6)                                      THEN

             JFAC(1) = 2
             JFAC(2) = 3

         ELSE IF (N.EQ.4)                                      THEN

             JFAC(1) = 2
             JFAC(2) = 2

         ELSE

*            * ABORT. THERE SHOULD BE AT LEAST TWO FACTORS.

             IFAC(1) = -1

             RETURN

         END IF

         NU = 1
         K  = 2

         GOTO 030

      END IF

***    FIRST, FIND PRIME FACTORS FACTORS OF N. 
  
      NU   = N
      K    = 0 
      L    = 1 
      IFAX = LFAC(L)

  010 CONTINUE

          IF (MOD(NU,IFAX).EQ.0)                               THEN 
              K = K+1 
              IF (K.GE.MAXFAC)                                 THEN
                  WRITE(*,1002)
                  GOTO 999
              END IF
              JFAC(K) = IFAX
              NU      = NU/IFAX
              IF (NU.EQ.1) GOTO 030 
                           GOTO 010
          END IF

  020     IFAX = LFAC(L+1)
          L    = L+1 

      IF (IFAX.GT.1) GOTO 010 

***    NOW CONSIDER NON PRIME FACTORS.

      CALL NOPRIM( NU,JFAC,K,MAXFAC )

  999 IF (NU.NE.1 .OR. K.GE.MAXFAC)                            THEN

***        ILLEGAL FACTOR OR TOO MANY FACTORS.

          IFAC(1) = -1
          RETURN

       END IF

*-----------------------------------------------------------------
***    NOW SORT FACTORS FROM SMALLEST ONWARD
***    AND STORE NF AS IFAC(1).
  
  030 CONTINUE

      IER = 0
      NF  = K

      DO 040 I=1,NF-1
         DO 040 L=I+1,NF
            IF (JFAC(I).GT.JFAC(L))                            THEN
                HOLD    = JFAC(L)
                JFAC(L) = JFAC(I)
                JFAC(I) = HOLD
            END IF
  040 CONTINUE

      IFAC(1) = NF
      DO 050 I=1,NF
         IFAC(I+1) = JFAC(I)
  050 CONTINUE

      WRITE(6,1001) (IFAC(I),I=1,NF+1)

      RETURN

*-------------------------------------------------------------------
 1001 FORMAT(' STOCKHAM  NFAX=',I4,', factors ...',19I4)
 1002 FORMAT(' Found too many factors... ',I4)

      END

* *** NOPRIM FACTORS WITH NON-PRIME FACTORS FOR SFFT.

      SUBROUTINE NOPRIM (N,JFAC,K,MAXFAC)

      IMPLICIT   none

***    AUTEUR:  B.Dugas - 13 novembre 1991.

      INTEGER    NBRNON,I,J,L,NU,IFAX,SETN
      PARAMETER( NBRNON = 10 )

      INTEGER    N,MAXFAC,K,K0(NBRNON),JFAC(MAXFAC),LFAC(NBRNON+1),
     +           TFAC(100,NBRNON)

      REAL*8     VALUE(NBRNON),NOVAL,SUM,VALM
      PARAMETER( NOVAL = 9999.99 )

      DATA       LFAC / 16,12,10,9,8,6,5,4,3,2,1 /
*-------------------------------------------------------------------
#     if defined (NEC)
*vdir novector
#     endif
      DO  I=1,NBRNON

          VALUE(I) = NOVAL

***        FIND FISRT PRIME FACTORS FACTORS OF N. 

          NU    = N
          K0(I) = 0
          L     = I
          IFAX  = LFAC(L)

  010     CONTINUE

              IF (NU.LE.12 .OR. NU.EQ.16)                      THEN

                  K0(I)         = K0(I)+1
                  TFAC(K0(I),I) = NU
                  NU            = 1

                  GOTO 020

              ELSE IF (MOD(NU,IFAX).EQ.0)                      THEN 

                  K0(I)         = K0(I)+1 
                  IF (K0(I)+K.GE.MAXFAC) CYCLE
                  TFAC(K0(I),I) = IFAX
                  NU            = NU/IFAX

                  IF (NU.EQ.1) GOTO 020 
                               GOTO 010

              END IF

              IFAX = LFAC(L+1)
              L    = L+1 

          IF (IFAX.GT.1) GOTO 010 

***        NU SHOULD BE ONE AT THIS POINT. 
***        STOP THIS PASS IF IT IS NOT.

  020     IF (NU.NE.1) CYCLE

***        CALCULATE VALUE FOR THIS SET.

          SUM = 0.0

          DO J=1,K0(I)
              SUM = SUM+1.0/FLOAT( TFAC(J,I) )
          END DO

          VALUE(I) = SUM

      END DO

***    CHOSE MINIMUM VALUE SET.

      VALM = VALUE(1)
      SETN = 1
      DO J=2,NBRNON
         IF (VALM.GT.VALUE(J))                                 THEN
             VALM = VALUE(J)
             SETN = J
         END IF
      END DO

      IF (VALM.NE.NOVAL)                                       THEN

***        TRANSFER INFO INTO JFAC,K AND UPDATE N.

          N = 1
          DO J=1,K0(SETN)
             JFAC(K+J) = TFAC(J,SETN)
          END DO
          K = K+K0(SETN)
         
      END IF

      RETURN
      END
      BLOCK DATA DATA_STKFFT8

!     Initialise variables of the COMMONs that need to be.

      INTEGER    MAXFAC
      PARAMETER( MAXFAC = nombre_de_facteurs )

      CHARACTER  STRING*8
      INTEGER    N, IFAC
      REAL*8,    POINTER :: TRIG(:,:)

      COMMON   / FTSCAL  / STRING
      COMMON   / STKFFT8 / TRIG, IFAC(MAXFAC+1), N

      DATA       STRING /  ' '   /
      DATA       N      /  -1    /

      END BLOCK DATA DATA_STKFFT8
