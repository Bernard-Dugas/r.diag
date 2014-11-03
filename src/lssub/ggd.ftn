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
      SUBROUTINE GGDX2(DX,GG,NLON,NLAT,RAD,CL)

      IMPLICIT    none

***    AUTHOR: NOV 13/80 - J.D.HENDERSON

***    SETS GRID DX TO THE LONGITUDE DERIVATIVE OF GRID GG.
***       I.E.   1./COS(LAT) * D(GG)/D(LON)

***    OTHER INPUT PARAMETERS...
***     RAD  CONTAINS LATITUDES IN RADIANS OF THE ROWS OF GG.
***     CL   CONTAINS THE CORRESPONDING COSINES.
***     NLAT NUMBER OF LATITUDES.
***     NLON NUMBER OF LONGITUDES.

      INTEGER     MAXJ
      PARAMETER ( MAXJ = 999 )

      INTEGER     I,J,NLON,NLONM1,NLONM2,NLAT
      REAL*8      RAD(NLAT),CL(NLAT),CLLIM,PI
      REAL        DX(NLON,NLAT),GG(NLON,NLAT),
     +            DELD,DELR,DELCOS

*--------------------------------------------------------------------
                              NLONM1 = NLON-1
      IF (MOD( NLON,2 ).EQ.0) NLONM1 = NLON
                              NLONM2 = NLONM1-1

      PI    = 4.0 * ATAN( 1.D0 )
      DELD  = 360./FLOAT( NLONM1 )
      DELR  = 2.*DELD*( PI/180. )

***    THE LIMIT ON CL IS DEFINED WITH
***    RESPECT TO A GRID WITH MAXJ LATITUDES.

      CLLIM = SIN( PI / MAXJ / 4.0 )

***    COMPUTE THE LONGITUDE DERIVATIVE.
***    USE CENTERED DIFFERENCES EVERYWHERE.

      DO 200 J=1,NLAT

***        CHECK FOR LATITUDE VERY CLOSE TO THE POLE.

          IF (CL(J).GE.CLLIM)                                  THEN

              DELCOS       = 1./( CL(J)*DELR )

              DO  I=2,NLONM2
                  DX(I ,J) = ( GG(I+1,J)-GG(I-1,J) )
     +                     *         DELCOS
              END DO

              DX(  1   ,J) = ( GG(2,J)-GG(NLONM1,J) )
     +                     *         DELCOS
              DX(NLONM1,J) = ( GG(1,J)-GG(NLONM2,J) )
     +                     *         DELCOS

          ELSE

              DO  I=1,NLON
                  DX(I, J) = 0.0
              END DO

          END IF

          IF (NLON.NE.NLONM1)
     +    DX(NLON  ,J)  =   DX(1,J)

  200 CONTINUE

      RETURN
      END

      SUBROUTINE GGDX3(DX,GG,NLON,NLAT,RAD,CL,ALON )

      IMPLICIT    none

      INTEGER     NLON,NLAT
      REAL*8      RAD(NLAT),CL(NLAT)
      REAL        DX(NLON,NLAT),GG(NLON,NLAT),ALON(NLON)

***    AUTHOR: SEP 11/03 - B. DUGAS (BASED ON GGDX2)

***    SETS GRID DX TO THE LONGITUDE DERIVATIVE OF GRID GG.
***       I.E.   1./COS(LAT) * D(GG)/D(LON)

***    OTHER INPUT PARAMETERS...
***     RAD  CONTAINS LATITUDES IN RADIANS OF THE ROWS OF GG.
***     CL   CONTAINS THE CORRESPONDING COSINES.
***     ALON CONTAINNS THE LONNGITUDES OF THE COLUMNS OF GG.
***     NLAT NUMBER OF LATITUDES.
***     NLON NUMBER OF LONGITUDES.

      INTEGER     I,J,NLM1
      LOGICAL     CYCLE,REPEAT
      REAL        DELTA,DEUXPI,OVCL

*--------------------------------------------------------------------

      DEUXPI = 8.0 * ATAN( 1.0 )
      NLM1   = NLON-1

***    CHECK FOR PERIODICITY IN LONGITUDE.

      DELTA = ( ALON(NLON) - ALON(NLM1) + ALON(2) - ALON(1) ) / 2.0

      REPEAT = .FALSE.
      CYCLE  = .FALSE.

      IF
     +(MOD( ABS( ALON(NLON)-ALON(1))      ,DEUXPI ).LT. 0.1*DELTA ) THEN
***        ALON(1) AND ALON(NLON) ARE "CLOSE".
          REPEAT = .TRUE.
      ELSE IF
     +(MOD( ABS( ALON(NLON)-ALON(1)+DELTA),DEUXPI ).LT. 0.1*DELTA ) THEN
***        ALON(1) AND ALON(NLON)+DELTA ARE "CLOSE".
          CYCLE = .TRUE.
      END IF

      DELTA = DELTA * 2.0

***    COMPUTE THE LONGITUDE DERIVATIVE.
***    USE CENTERED DIFFERENCES NEARLY EVERYWHERE.

      DO  J=1,NLAT

          OVCL = 1./CL(J)

          DO  I=2,NLM1
              DX(  I ,J) = ( GG(I+1,J) - GG(I-1,J) )   * OVCL
     +                   / ( ALON(I+1) - ALON(I-1) )
          END DO

          IF (REPEAT)                                          THEN
              DX(  1 ,J) = ( GG(  2,J) - GG(NLM1,J) )  * OVCL
     +                   /           DELTA
              DX(NLON,J) = DX(1,J)
          ELSE IF (CYCLE)                                      THEN
              DX(  1 ,J) = ( GG(  2,J) - GG(NLON,J) )  * OVCL
     +                   /           DELTA
              DX(NLON,J) = ( GG(  1,J) - GG(NLM1,J) )  * OVCL
     +                   /           DELTA
          ELSE
              DX(  1 ,J) = ( GG(  2,J) - GG(   1,J) )  * OVCL
     +                   /   ( ALON(2) - ALON(1) )
              DX(NLON,J) = ( GG(NLON,J) - GG(NLM1,J) ) * OVCL
     +                   / ( ALON(NLON) - ALON(NLM1) )
          END IF

      END DO

      RETURN
      END

      SUBROUTINE GGDY2(DY,GG,NLON,NLAT,RAD)

***    AUTHOR: NOV 13/80 - J.D.HENDERSON

***    SETS GRID DY TO THE  LATITUDE DERIVATIVE OF GRID GG.

***    OTHER INPUT PARAMETERS...
***     RAD  CONTAINS LATITUDES IN RADIANS OF THE ROWS OF GG.
***     NLAT NUMBER OF LATITUDES.
***     NLON NUMBER OF LONGITUDES.

      IMPLICIT none

      INTEGER  I,J, NLAT,NLATM,NLON
      REAL     DY(NLON,NLAT),GG(NLON,NLAT)
      REAL*8   RAD(NLAT)

*--------------------------------------------------------------------
      NLATM = NLAT-1

***    COMPUTE THE  LATITUDE DERIVATIVE.
***    ENDS USE ONE-SIDED DIFFERENCES. INTERIOR USES CENTERED DIFF.

      DO 100 J=2,NLATM
          DO 100 I=1,NLON

              DY(I,J) = ( GG(I,J+1)-GG(I,J-1) )
     +                /  ( RAD(J+1)-RAD(J-1) )

  100 CONTINUE

      DO 200 I=1,NLON

          DY(I,1)    =    ( GG(I,2)-GG(I,1) )
     +               /     ( RAD(2)-RAD(1) )

          DY(I,NLAT) = ( GG(I,NLAT)-GG(I,NLATM) )
     +               /  ( RAD(NLAT)-RAD(NLATM) )
  200 CONTINUE

      RETURN
      END
