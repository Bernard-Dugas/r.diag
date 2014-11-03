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
      SUBROUTINE HTEST (P,NWDS,ALPHA,IKIND,K,D)

      IMPLICIT  REAL*8 (A-H,O-Z), INTEGER (I-N)
     
***    F. MAJAESS -- SEPT. 30, 1987 (IN COMPUTING THE FIELDS OF
***                                 TRANSFORMED SIGNIFICANCE LEVELS K, 
***                                 USE LOG BASE 5 INSTEAD OF LOG BASE 
***                                 2. I.E USE P=ALPHA/5**(K-1) INSTEAD
***                                 OF P=ALPHA/2**(K-1)).
***    F. ZWIERS  -- NOV.  27, 1984

***    THIS SUBROUTINE CONDUCTS A TEST OF HYPOTHESIS BY EXAMINING P-VALUES. 

***    INPUT:  
***       P    : ARRAY OF P-VALUES.
***       NWDS : LENGTH OF ARRAY P.
***       ALPHA: THE SIGNIFICANCE LEVEL OF THE TEST. 
***       IKIND=-1 REJECT THE NULL HYPOTHESIS IF P > 1 - ALPHA 
***                (IE - CONDUCT A ONE TAILED TEST AND REJECT IF 
***                      THE TEST STATISTIC IS UNUSUALLY SMALL). 
***            = 0 INDICATES A TWO TAILED TEST (REJECT IF
***                P > 1 - ALPHA/2   OR IF  P < ALPHA/2) 
***                (IE - REJECT IF THE TEST STATISTIC IS 
***                      UNUSUALLY SMALL OR LARGE).
***            = 1 REJECT THE NULL HYPOTHESIS IF P < ALPHA 
***                (IE - CONDUCT A ONE TAILED TEST AND REJECT IF 
***                      THE TEST STATISTIC IS UNUSUALLY LARGE). 

***    OUTPUT: 
***       K    : REAL ARRAY OF TRANSFORMED SIGNIFICANCE LEVELS.
***              THE K-VALUES ARE LOG TO BASE 5 OF ALPHA DEVIDED BY
***              THE SIGNIFICANCE OF THE OBSERVED STATISTIC. THUS A
***              K-VALUE OF 1 INDICATES THAT THE OBSERVED STATISTIC
***              IS JUST SIGNIFICANT AT THE ALPHA SIGNIFICANCE LEVEL.
***              IN GENERAL K-VALUE INDICATES THAT THE OBSERVED
***              STATISTIC WOULD BE SIGNIFICANT AT THE 
***                             ALPHA/5**(K-1) 
***              SIGNFICANCE LEVEL.
***              IF PLOTTED WITH UNIT CONTOUR INTERVALS, SUCCESSIVE
***              CONTOURS WILL ENCLOSE REGIONS WHERE LOCALLY IT IS 
***              FIVE TIMES AS UNLIKELY THAT VALUES OF THE OBSERVED
***              STATISTICS ARE CONSISTENT WITH THE NULL HYPOTHESIS
***              THAN IN REGIONS OUTSIDE THE NEXT LOWER CONTOUR. 
***       D    : ARRAY CONSISTING OF 0'S (HYPOTHESIS ACCEPTED) 
***              AND 1'S (HYPOTHESIS REJECTED).

      INTEGER   IKIND,NWDS
      DIMENSION P(NWDS),K(NWDS),D(NWDS)

      REAL*8    LN5,LNABS5,K

      DATA      EPS / 1.E-12 /

C---------------------------------------------------------------------------- 
      LN5    = LOG( 5.0D0 ) 
      LNABS5 = LOG( ALPHA )/LN5

***    CONVERT THE P-VALUES TO SIGNFICANCE LEVELS AND STORE THEM 
***    IN ARRAY D WHICH IS USED TEMPORARILY AS A WORKING ARRAY.

      IF (IKIND.EQ.-1)                                         THEN 
          DO  I=1,NWDS 
              D(I) = 1.0-P(I) 
          END DO
      ELSE IF (IKIND.EQ.0)                                     THEN
          DO  I=1,NWDS 
              D(I) = 2.0*MIN( P(I),1.0-P(I) ) 
          END DO
      ELSE
          DO  I=1,NWDS 
              D(I) = P(I) 
          END DO
      END IF 

***    BEFORE COMPUTING THE K-VALUES, MAKE SURE THE SIGNIFICANCE LEVEL 
***    IS NOT ZERO. (BECAUSE WE PLAN TO TAKE LOGS).

      DO  I=1,NWDS
          D(I) = MAX( D(I),EPS ) 
      END DO

***    COMPUTE THE TRANSFORMED SIGNIFICANCE LEVELS.(K-VALUES). 

      DO  I=1,NWDS
          K(I) = LNABS5-LOG( D(I) )/LN5+1.0 
      END DO

***    MAKE DECISION ABOUT WHETHER OR NOT TO REJECT. 
***    (I.E. COMPUTE D-VALUES).

      DO  I=1,NWDS
          D(I) = 0.0 
      END DO
      DO  I=1,NWDS
          IF (K(I).GE.1.0)                                     THEN
              D(I) = 1.0
          ELSE
              D(I) = 0.0
          END IF
      END DO

      RETURN
C---------------------------------------------------------------------------- 

      END 
      SUBROUTINE PROBF (F,V1,V2,NWDS,P,THETA,Z,B1) 

      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
     
***    F. ZWIERS  -- JAN 18/85 

***    EXACT CALCULATION OF PROB( F > EF ) WHERE A HAS V1 AND V2
***    DEGREES OF FREEDOM.

***    ARGUEMENT LIST: 
***     INPUT -
***       F    - ARRAY OF F VALUES 

***              ***************************************************** 
***              *                                                   * 
***              *  THE CONTENTS OF ARRAY F ARE MODIFIED BY THIS     * 
***              *  SUBROUTINE                                       * 
***              *                                                   * 
***              ***************************************************** 

***       V1   - DEGREES OF FREEDOM FOR THE NUMERATOR OF THE
***              F-RATIOS (SCALAR)
***       V2   - DEGREES OF FREEDOM FOR THE DENOMINATOR OF THE
***              F-RATIOS (SCALAR) 
***       NWDS - LENGTH OF ARRAY F 

***     OUTPUT - 
***       P    - ARRAY OF PROBABILITIES OF OBSERVING F-RATIOS WITH V1
***              AND V2 DEGREES OF FREEDOM WHICH ARE LARGER THAN THE 
***              CORRESPONDING F-VALUES CONTAINED IN ARRAY F.

***     WORK SPACE - 
***       THETA,Z,B1 - WORK ARRAYS OF LENGTH NWDS

***    REF - LACKRITZ, J.R., 1984: EXACT P VALUES FOR F AND T TESTS, 
***          THE AMERICAN STATISTICIAN, 38, 312-314. 

      INTEGER   V1,V2,NWDS
      DIMENSION F(NWDS),P(NWDS),THETA(NWDS),Z(NWDS),B1(NWDS)
      LOGICAL   N1EVEN,N2EVEN 
      DATA      EPS/1.E-6/ 

C-----------------------------------------------------------------------------
      PI = 4.0*ATAN( 1.D0 )

      IF (V1.LT.1 .OR. V2.LT.1)                                THEN
          WRITE(6,6010) V1,V2
          CALL                                     XIT('  Probf ',-1 ) 
      ENDIF 

      DO  I=1,NWDS 
         Z(I) = MAX( EPS,F(I) ) 
         F(I) = Z(I)
      END DO

      N1 = V1 
      N2 = V2 
      IF (V2.LT.V1)                                            THEN
          DO  I=1,NWDS
              F(I) = 1.0/F(I)
          END DO
          N1 = V2
          N2 = V1
      END IF 

      FN1 = FLOAT( N1 ) 
      FN2 = FLOAT( N2 ) 
      M1  = N1/2 
      M2  = N2/2 

                      N1EVEN =.TRUE.
      IF (M1*2.LT.N1) N1EVEN =.FALSE. 
                      N2EVEN =.TRUE.
      IF (M2*2.LT.N2) N2EVEN =.FALSE. 

      FN1BN2 = FN1/FN2
      DO  I=1,NWDS 
          THETA(I)=F(I)*FN1BN2
      END DO

      IF (N1EVEN)                                              THEN

***    V1 EVEN 

          DO  I=1,NWDS
              F(I) = THETA(I)/(1.0 + THETA(I))
              Z(I) = 1.0 
              P(I) = 1.0 
          END DO
          IF (M1.GT.1)                                         THEN
              DO  I=1,M1-1 
                  FTWOI = FLOAT( I+I ) 
                  FACTI = (FN2 + FTWOI - 2.0)/FTWOI
                  DO  J=1,NWDS
                      Z(J) = Z(J)*F(J)*FACTI
                      P(J) = P(J) + Z(J)
                  END DO 
              END DO
          END IF
          POWER=FN2/2.0
          DO  I=1,NWDS
              P(I)=P(I)/((1.0 + THETA(I))**POWER) 
          END DO

      ELSE IF (N2EVEN)                                         THEN

***    V1 ODD, V2 EVEN 

          DO  I=1,NWDS
              THETA(I) = 1.0/THETA(I) 
              Z(I)     = 1.0
              P(I)     = 1.0
              F(I)     = THETA(I)/(1.0 + THETA(I))
          END DO

          IF (M2.GT.1)                                         THEN
              DO  I=1,M2-1 
                  FTWOI = FLOAT( I+I ) 
                  FACTI = (FN1 + FTWOI - 2.0)/FTWOI
                  DO  J=1,NWDS
                      Z(J) = Z(J)*F(J)*FACTI
                      P(J) = P(J) + Z(J)
                  END DO
              END DO
          END IF

          POWER = FN1/2.0
          DO  I=1,NWDS
              P(I) = 1.0 - P(I)/((1.0 + THETA(I))**POWER) 
          END DO

      ELSE

***    V1 ODD, V2 ODD

          IF (N1.EQ.1 .AND. N2.EQ.1)                           THEN

              DO  I=1,NWDS 
                  P(I) = 2.0*ATAN( 1.0/SQRT( THETA(I) ) )/PI 
              END DO

          ELSE 

              DO  I=1,NWDS 
                  F(I)=2.0/(1.0 + THETA(I))
                  Z(I)=F(I)
                  P(I)=F(I)
              END DO

              IF (M2.GT.1)                                     THEN 
                  DO  I=2,M2
                      FTWOI = FLOAT( I+I )
                      FACTI = FLOAT( I-1 )/(FTWOI-1.0)
                      DO  J=1,NWDS 
                          Z(J) = Z(J)*F(J)*FACTI 
                          P(J) = P(J) + Z(J) 
                      END DO
                  END DO
              END IF 

              DO  I=1,NWDS 
                  B1(I)=SQRT(THETA(I))*P(I)/PI 
              END DO

              IF (N1.EQ.1)                                     THEN 
                  DO  I=1,NWDS
                      P(I) = 2.0*ATAN( 1.0/SQRT( THETA(I) ) )/PI - B1(I)
                  END DO
              ELSE
                  FACT = 1.0 
                  DO  I=1,M2-1
                      FACT = FACT*FLOAT( I )/(FLOAT( I )+0.5) 
                  END DO
                  FACT = 2.0*FACT
                  DO  I=1,NWDS
                      F(I) = 2.0*THETA(I)/(1.0 + THETA(I))
                      Z(I) = 1.0
                      P(I) = 0.0
                  END DO
                  DO  I=1,M1
                      FACTI = FLOAT( M2 + I - 1 )/FLOAT( I + I - 1 )
                      DO  J=1,NWDS 
                          Z(J) = Z(J)*F(J)*FACTI 
                          P(J) = P(J) + Z(J) 
                      END DO
                  END DO
                  DO  I=1,NWDS
                      P(I) =  P(I)*FACT 
     +                     / (PI*SQRT(THETA(I))*(1.0 + THETA(I))**M2) 
                      P(I) =  P(I) - B1(I) 
     +                     +  2.0*ATAN( 1.0/SQRT( THETA(I) ) )/PI 
                  END DO
              END IF 

          END IF

      END IF 

      IF (V2.LT.V1)                                            THEN 
          DO  I=1,NWDS
              P(I) = 1.0 - P(I) 
          END DO
      END IF

      RETURN
C-----------------------------------------------------------------------------

 6010 FORMAT(' Probf called with degrees of freedom ',2I5)

      END 
      SUBROUTINE PROBT (T,NDF,NWDS,P,WK1,WK2,WK3,WK4)

      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
     
***    F. ZWIERS  DEC 6,1984 

***    COMPUTE PROB( T > T-OBSERVED )

      INTEGER   NDF,NWDS
      DIMENSION T(NWDS),P(NWDS),WK1(NWDS),WK2(NWDS),WK3(NWDS),WK4(NWDS) 

C-----------------------------------------------------------------------------
      DO  I=1,NWDS
          WK1(I)=T(I)*T(I) 
      END DO

      CALL PROBF( WK1,1,NDF,NWDS,P,WK2,WK3,WK4 )

      DO  I=1,NWDS
          IF (T(I).GE.0.0)                                     THEN
              P(I) =    P(I)/2.0
          ELSE
              P(I) = 1.-P(I)/2.0
          END IF
      END DO

      RETURN
C-----------------------------------------------------------------------------

      END 
