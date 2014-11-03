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
C     $Log: interp.ftn,v $
C     Revision 3.15  2014/09/25 18:42:03  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.14  2006/04/26 16:12:29  dugas
C     Enlever les routines grid_to_grid, utiliser les versions dans rmnlib.
C
C     Revision 3.13  2005/07/28 17:24:12  dugas
C     Modifier le code pour enlever les messages d'avertissement de F90.
C
C     Revision 3.12  2005/04/12 16:45:06  dugas
C     Considerer le cas ILG=ILG1 dans IGGSL.
C
C     Revision 3.11  2004/11/22 03:44:35  dugas
C     Ne plus calculer SLON dans LLCAL.
C
C     Revision 3.10  2004/11/08 20:49:21  dugas
C     Ajouter les routine du groupe GRID_TO_GRID (provenant de MFV).
C
C     Revision 3.9  2002/04/22 13:31:18  dugas
C     Corriger le calcul de NLG dans GGILL2
C
C     Revision 3.8  1999/10/28 01:23:25  armnrbd
C     Verifier la concordance entre ILG et ILG1 avant de
C      sortir des routine LLEGG2 et LLIGG2.
C
C     Revision 3.7  1999/04/08 19:23:40  armnrbd
C     Enlever deux ENDIF en bout de ligne dans HRTOLR et HRALR.
C
C     Revision 3.6  1998/12/01 18:53:18  armnrbd
C     Corriger LLIGGM pour les coupes zonales.
C
C     Revision 3.5  1997/06/02  16:39:56  armnrbd
C     Ajouter la routine IGGSL.
C
C     Revision 3.5  1997/06/02  14:28:36  armnrbd
C     Ajouter la routine IGGSL.
C
C     Revision 3.4  1997/04/30  19:43:37  armnrbd
C     Ajouter le support des coupes zonales dans les routines LGINT*.
C
C     Revision 3.3  1995/11/23  02:38:45  armnrbd
C     Ajouter les routines GINTLM et LLIGGM.
C
C     Revision 3.2  1995/11/01  20:04:54  armnrbd
C     De-allouer les appels multi-taches.
C     Ajouter la routine GGIPS3.
C
C     Revision 3.1  1994/11/17  14:13:33  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:55:40  13:55:40  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.1  93/11/18  11:37:00  armnrbd
C     Corriger un bogue dans PSCAL (mauvaise declaration de D60)
C     
C     Revision 2.0  93/10/13  13:31:49  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.5  93/08/23  08:42:22  armnrbd
C     Corriger un enonce DOACCROSS.
C     
C     Revision 1.4  93/08/19  16:20:37  armnrbd
C     Modifications cosmetiques.
C     
C     Revision 1.3  93/07/08  14:06:32  armnrbd
C     Implanter version hemispherique de GGILL2/GGPNT2.
C     
C     Revision 1.2  92/11/26  17:47:04  armnrbd
C     BugFix.
C     
C     Revision 1.1  92/11/26  16:48:56  armnrbd
C     Ajouter les routines LLCAL, PSCALL, GGILL et XYFLL pour le
C        programme de calcul des sous-aires (SUBAREA).
C     
C     Revision 1.0  92/02/21  11:33:22  armnrbd
C     Initial revision
C     

      SUBROUTINE fmean2 (G,LI,LJ,AVRG,N) 

***   ****   JAN 1975  -  JOHN D. HENDERSON  **** 
***    CALCULATES MEAN OF G(LI,LJ) OMITTING N BORDER ROWS. 

      IMPLICIT none

      INTEGER  LI,LJ, IL,JL, IH,JH, I,J, N
      REAL     G(LI,LJ), PTS, AVRG
      REAL*8   SUM 

*----------------------------------------------------------------------- 
      IL  = 1+N 
      JL  = 1+N 
      IH  = LI-N 
      JH  = LJ-N 
      PTS = FLOAT( (IH-IL+1)*(JH-JL+1) )

      SUM = 0. 
      DO 20 J=JL,JH 
          DO 20 I=IL,IH 
              SUM  = SUM+G(I,J) 
   20 CONTINUE
      AVRG = SUM/PTS 

      RETURN 
*--------------------------------------------------------------------

      END 
      SUBROUTINE ggill (DAT,G, I1,I2,J1,J2,LX,LY) 
 
***    AUTHOR: AUG 10/87 - M.SUTCLIFFE, CCRN.

***    EXTRACT SEGMENTS OF ARRAY G AND PUT THEM INTO ARRAY DAT. 
***    I1,I2 AND J1,J2 ARE THE LEFT,RIGHT AND LOWER,UPPER GRID 
***    INDICES OF THE WINDOW. 

***    LX IN IS THE NUMBER OF INPUT GRID POINTS IN THE X DIRECTION. 
***    LX,LY OUT ARE THE NUMBERS OF GRID POINTS IN THE WINDOW. 
 
      IMPLICIT none

      INTEGER  I1,I2,J1,J2,LX,LY
      REAL     G(LX),DAT(LX) 

      INTEGER  I,J,IWINDO,IFULLG

*--------------------------------------------------------------------- 
***    EXTRACT WINDOW. 
 
      IWINDO = 0 
      DO 400 J=J1,J2 

          IF (I1.GE.I2)                                        THEN 

              DO 100 I=I1,LX-1 
                  IWINDO      = IWINDO+1 
                  IFULLG      = (J-1)*LX+I 
                  DAT(IWINDO) = G(IFULLG) 
  100         CONTINUE 
              DO 200 I=1,I2 
                  IWINDO      = IWINDO+1 
                  IFULLG      = (J-1)*LX+I 
                  DAT(IWINDO) = G(IFULLG) 
  200         CONTINUE 

          ELSE 

              DO 300 I=I1,I2 
                  IWINDO      = IWINDO+1 
                  IFULLG      = (J-1)*LX+I 
                  DAT(IWINDO) = G(IFULLG) 
  300         CONTINUE 

          END IF 

  400 CONTINUE 
 
      LY = J2-J1+1 
      IF (I1.GE.I2)                                            THEN 
          LX = I2-I1+LX 
      ELSE 
          LX = I2-I1+1 
      ENDIF 
 
      RETURN 
*--------------------------------------------------------------------

      END 
      SUBROUTINE ggill2 (LL,NLG1,NLAT,GG,ILG1,ILAT,SLAT, INTERP) 
 
***    MAY  1/80 - J.D.HENDERSON 

***    INTERPOLATES GLOBAL LAT-LONG GRID LL(NLG1,NLAT) 
***     FROM GLOBAL GAUSSIAN GRID GG(ILG1,ILAT). 

***    SLAT   = LAT (DEG) OF GAUSSIAN GRID ROWS FROM THE S POLE. 
***    INTERP = (1,3) FOR (LINEAR,CUBIC) INTERPOLATION. 
***             (OTHERWISE THE GRID GG IS SET TO ZERO). 
 
      IMPLICIT none

      INTEGER  NLG1,NLAT,ILG1,ILAT,INTERP
      REAL     LL(NLG1,NLAT),GG(ILG1,ILAT) 
      REAL     SLAT(ILAT) 

      INTEGER  NLG,I,J
      LOGICAL  GLOB,NHEM
      REAL     DX,DY,VAL,LAT0,
     +         DLAT,DLON,AVG

      SAVE     NP
      EXTERNAL mpserv
      INTEGER  mpserv,NP
CC$   INTEGER  Bidon

      EXTERNAL GGPNL2,GGPNT2,FMEAN2

*----------------------------------------------------------------------- 
***    CHECK FOR HEMISPHERIC OR GLOBAL INPUT DATA.

      IF (SLAT(ILAT).LT.90.)                                   THEN
          GLOB = .FALSE.
          NHEM = .FALSE.
          LAT0 =  00.
       ELSE IF (SLAT(1).GT.90.)                                THEN
          GLOB = .FALSE.
          NHEM = .TRUE.
          LAT0 =  90.
      ELSE
          GLOB = .TRUE.
          NHEM = .FALSE.
          LAT0 =  00.
      END IF

***    DEFINE CONSTANT INCREMENTS.

      NLG = NLG1-1

      DX  = 360./FLOAT( NLG ) 
      IF (GLOB)                                                THEN
          DY  = 180./FLOAT( NLAT-1 ) 
      ELSE
          DY  =  90./FLOAT( NLAT-1 ) 
      END IF
 
***    CHECK THE NUMBER OF PARALLEL THREADS AVAILABLE.

                   NP    = mpserv('THREADS',NP)
CC$     IF (NP.EQ.1) Bidon = mpserv('DESTROY',Bidon)

***    DO THE INTERPOLATION

CC$doacross local(I,J,DLAT,DLON,VAL)

      DO 100 J=1,NLAT 
          DLAT=LAT0+FLOAT(J-1)*DY 
          DO I=1,NLG 
              DLON = FLOAT(I-1)*DX 
              VAL      = 0.0
              IF (INTERP.EQ.1) 
     +            CALL GGPNL2( VAL, GG,ILG1,ILAT,DLAT,DLON,SLAT ) 
              IF (INTERP.EQ.3) 
     +            CALL GGPNT2( VAL, GG,ILG1,ILAT,DLAT,DLON,SLAT ) 
              LL(I,J) = VAL 
          END DO
          LL(NLG1,J)  = LL(1,J) 
  100 CONTINUE
 
CC$    Bidon = mpserv('BLOCK',Bidon)

***    AS NEED BE...
 
      IF (GLOB .OR. NHEM)                                      THEN

***        SET N POLE TO MEAN OF TOP ROW. 

          CALL FMEAN2( LL(1,NLAT),NLG,1,AVG, 0 ) 
          DO  I=1,NLG1 
              LL(I,NLAT) = AVG 
          END DO

      END IF
      IF (GLOB .OR. .NOT.NHEM)                                 THEN

***        SET S POLE TO MEAN OF BOTTOM ROW. 

          CALL FMEAN2( LL(1,1),NLG,1,AVG, 0 ) 
          DO  I=1,NLG1 
              LL(I,1) = AVG 
          END DO

      END IF
 
      RETURN 
*--------------------------------------------------------------------

      END 
      SUBROUTINE ggips2 (PS,LI,LJ,IP,JP,D60,DGRW,NHEM, 
     +                   GG,ILG1,ILAT,ANG,INTERP) 

***    MAI 15/91 - B.Dugas, RPN. (Version parrallele SGI)
***    DEC  1/80 - J.D.HENDERSON 

***    INTERPOLATES POLAR STEREOGRAPHIC GRID PS(LI,LJ) 
***    FROM GLOBAL GAUSSIAN GRID GG(ILG1,ILAT). 
***    ANG CONTAINS DEGREES LATITUDE FROM THE SOUTH POLE. 
***    GGPNL2 AND GGPNT2 BOTH REQUIRE 
***    LATITUDE BETWEEN 0 AND 180 AND LONGITUDE BETWEEN 0 AND 360. 

      IMPLICIT none

      INTEGER  LI,LJ,IP,JP,NHEM,ILG1,ILAT,INTERP, I,J
      REAL     PS(LI,LJ),GG(ILG1,ILAT), ANG(ILAT), D60,DGRW, 
     +         X,Y, DLAT,DLON, VAL

      SAVE     NP
      EXTERNAL mpserv
      INTEGER  mpserv,NP
CC$   INTEGER  Bidon

      EXTERNAL LLFXY,GGPNL2,GGPNT2

*-------------------------------------------------------------------- 
***    CHECK THE NUMBER OF PARALLEL THREADS AVAILABLE.

                   NP    = mpserv('THREADS',NP)
CC$    IF (NP.EQ.1) Bidon = mpserv('DESTROY',Bidon)

CC$doacross local(I,J,DLAT,DLON,X,Y,VAL)
 
      DO 320 J=1,LJ 
          Y = FLOAT( J-JP ) 
          DO 310 I=1,LI 
              X    = FLOAT( I-IP ) 
              CALL LLFXY( DLAT,DLON, X,Y,D60,DGRW,NHEM ) 

                              DLAT = DLAT+90. 
              IF (DLON.LT.0.) DLON = 360.+DLON 

              VAL  = 0. 
              IF (INTERP.EQ.1)
     +            CALL GGPNL2( VAL, GG,ILG1,ILAT,DLAT,DLON,ANG ) 
              IF (INTERP.EQ.3) 
     +            CALL GGPNT2( VAL, GG,ILG1,ILAT,DLAT,DLON,ANG ) 
              PS(I,J) = VAL 
  310     CONTINUE
  320 CONTINUE
 
CC$    Bidon = mpserv('BLOCK',Bidon)

      RETURN 
*--------------------------------------------------------------------

      END 
      SUBROUTINE ggips3 (PS,LI,LJ,PI,PJ,D60,DGRW,NHEM, 
     +                   GG,ILG1,ILAT,ANG,INTERP) 

***    MAI 15/91 - B.Dugas, RPN. (Version parrallele SGI)
***    DEC  1/80 - J.D.HENDERSON 

***    INTERPOLATES POLAR STEREOGRAPHIC GRID PS(LI,LJ) 
***    FROM GLOBAL GAUSSIAN GRID GG(ILG1,ILAT). 
***    ANG CONTAINS DEGREES LATITUDE FROM THE SOUTH POLE. 
***    GGPNL2 AND GGPNT2 BOTH REQUIRE 
***    LATITUDE BETWEEN 0 AND 180 AND LONGITUDE BETWEEN 0 AND 360. 

      IMPLICIT none

      INTEGER  LI,LJ,NHEM,ILG1,ILAT,INTERP, I,J
      REAL     PS(LI,LJ),GG(ILG1,ILAT), ANG(ILAT), 
     +         PI,PJ,D60,DGRW, X,Y, DLAT,DLON, VAL

      SAVE     NP
      EXTERNAL mpserv
      INTEGER  mpserv,NP
CC$   INTEGER  Bidon

      EXTERNAL LLFXY,GGPNL2,GGPNT2

*-------------------------------------------------------------------- 
***    CHECK THE NUMBER OF PARALLEL THREADS AVAILABLE.

                   NP    = mpserv('THREADS',NP)
CC$    IF (NP.EQ.1) Bidon = mpserv('DESTROY',Bidon)

CC$doacross local(I,J,DLAT,DLON,X,Y,VAL)
 
      DO 320 J=1,LJ 
          Y = J-PJ
          DO 310 I=1,LI 
              X    = I-PI
              CALL LLFXY( DLAT,DLON, X,Y,D60,DGRW,NHEM ) 

                              DLAT = DLAT+90. 
              IF (DLON.LT.0.) DLON = 360.+DLON 

              VAL  = 0. 
              IF (INTERP.EQ.1)
     +            CALL GGPNL2( VAL, GG,ILG1,ILAT,DLAT,DLON,ANG ) 
              IF (INTERP.EQ.3) 
     +            CALL GGPNT2( VAL, GG,ILG1,ILAT,DLAT,DLON,ANG ) 
              PS(I,J) = VAL 
  310     CONTINUE
  320 CONTINUE
 
CC$    Bidon = mpserv('BLOCK',Bidon)

      RETURN 
*--------------------------------------------------------------------

      END 
      SUBROUTINE ggpnt2 (VAL,GG,ILG1,ILAT,DLAT,DLON,WRKS) 

***    APR 21/83 - R.LAPRISE. 
***     WRITTEN BY ROGER DALEY - NOVEMBER 1975. 

***     THIS SUBROUTINE INTERPOLATEDS FROM A GLOBAL GAUSSIAN GRID 
***     GG(ILG1,ILAT) TO A SINGLE POINT WITH LATITUDE (DEGREES) DLAT 
***     AND LONGITUDE DLON. DLAT IS MEASURED FROM THE SOUTH POLE 
***     AND DLON IS MEASURED EASTWARD FROM GREENWICH AND MUST SATISFY 
***     (0.LE.DLON.LE.360).  BECAUSE THE GAUSSIAN LATITUDES ARE ONLY
***     APPROXIMATELY EQUALLY SPACED, SOME SEARCHING WILL BE DONE IN
***     THE GENERAL CASE. WRKS CONTAINS THE LATITUDES IN DEGREES FROM
***     THE S. POLE. 

      IMPLICIT  none

      INTEGER   ILG,ILG1,ILAT,IRX,N
      REAL      VAL, GG(ILG1,ILAT), DLAT,DLON,      DX,DY,RX,RY, 
     +          RS1,RS2,RSOUTH,     RN1,RN2,RNORTH, LATI,LATF,
     +          AG,AGM,AGP,         WRKS(2) 

      EXTERNAL  GINTC2,FMEAN2

*----------------------------------------------------------------------- 
***    DEFINE LATITUDINAL BOUNDARIES.

      IF (WRKS(ILAT).LT.90.)                                   THEN
          LATF =  090.
          LATI =  000.
       ELSE IF (WRKS(1).GT.90.)                                THEN
          LATF =  180.
          LATI =  090.
      ELSE
          LATF =  180.
          LATI =  000.
      END IF

***    DEFINE DX,DY (DELTAS) AND RX,RY (POSITIONS).

                               ILG = ILG1-1
      IF (MOD( ILG1,2 ).EQ.0 ) ILG = ILG1

      DX  = FLOAT( ILG )/360. 
      RX  = DLON*DX + 1. 
      IRX = INT( RX ) 

      DY  = 1./(WRKS(2)-WRKS(1)) 
      IF (DLAT.GT.(LATI+LATF)/2.)                              THEN
          RY = (DLAT-WRKS(ILAT))*DY + FLOAT( ILAT )
      ELSE
          RY = (DLAT-WRKS( 01 ))*DY + 1. 
      END IF

      IF (DLAT.GT.WRKS(000002) .AND.
     +    DLAT.LT.WRKS(ILAT-1) )                               THEN
 
              N   = INT( RY ) 
              AG  = WRKS(N) 
              AGP = WRKS(N+1) 
              AGM = WRKS(N-1) 

          IF (DLAT.LT.AG)                                      THEN
              RY  = FLOAT( N-1 ) + (DLAT-AGM)/(AG -AGM) 
          ELSE
              RY  = FLOAT(  N  ) + (DLAT-AG )/(AGP-AG ) 
          END IF

      END IF

      IF (DLAT.GT.WRKS(1)   .AND.
     +    DLAT.LT.WRKS(ILAT))                                  THEN

***        DO THE INTERPOLATION IN THE GENERAL CASE.

          CALL GINTC2( VAL, GG,ILG1,ILAT, RX,RY ) 

      ELSE IF (DLAT.LE.WRKS(1))                                THEN

***        SPECIAL CASE FOR DLAT SOUTH OF FIRST GAUSSIAN LATITUDE. 

          CALL FMEAN2( GG(1,1),ILG1,1,RS1,0 ) 
          CALL FMEAN2( GG(1,2),ILG1,1,RS2,0 ) 
          RSOUTH =  RS1 + (RS2-RS1)*DY*(LATI-WRKS(1))
          VAL    = (RX-FLOAT( IRX   ))*GG(IRX  ,1)
     +           - (RX-FLOAT( IRX+1 ))*GG(IRX+1,1) 
          VAL    = ( (WRKS(0001)-DLAT)*RSOUTH + (DLAT-LATI)*VAL )
     +           /           ( WRKS(0001)-LATI )

      ELSE IF (DLAT.GE.WRKS(ILAT))                             THEN

***        SPECIAL CASE FOR DLAT NORTH OF LAST GAUSSIAN LATITUDE. 

          CALL FMEAN2( GG(1,ILAT),  ILG1,1,RN1,0 ) 
          CALL FMEAN2( GG(1,ILAT-1),ILG1,1,RN2,0 ) 
          RNORTH =  RN1 + (RN1-RN2)*DY*(LATF-WRKS(ILAT))
          VAL    = (RX-FLOAT( IRX   ))*GG(IRX  ,ILAT)
     +           - (RX-FLOAT( IRX+1 ))*GG(IRX+1,ILAT) 
          VAL    = ( (DLAT-WRKS(ILAT))*RNORTH + (LATF-DLAT)*VAL )
     +           /           ( LATF-WRKS(ILAT) ) 

      END IF

      RETURN 
*--------------------------------------------------------------------

      END 
      SUBROUTINE ggpnl2 (VAL,GG,ILG1,ILAT,DLAT,DLON,SLAT) 

***    MAY  1/80 - J.D.HENDERSON 

***    LINEAR INTERPOLATION AT POINT (DLON,DLAT) 
***    IN GLOBAL GAUSSIAN GRID GG(ILG1,ILAT). 

      IMPLICIT none

      INTEGER  ILG,ILG1,ILAT,I,J,K
      REAL     VAL, BOT,TOP, X,Y, DI,FI,
     +         GG(ILG1,ILAT), SLAT(ILAT),
     +         DLAT,DLON

*----------------------------------------------------------------------- 
                               ILG = ILG1-1
      IF (MOD( ILG1,2 ).EQ.0 ) ILG = ILG1

      DI = 360./FLOAT( ILG ) 
      FI = DLON/DI+1. 
      I  = INT( FI ) 
      X  = FI-FLOAT( I ) 

      DO 110 K=1,ILAT 
          J=K-1 
          IF (DLAT.LT.SLAT(K)) GOTO 120 
  110 CONTINUE

  120 IF (J.EQ.0) J=1 

      IF (DLAT.GT.SLAT(ILAT)) J=ILAT-1 
      Y   = (DLAT-SLAT(J)) / (SLAT(J+1)-SLAT(J)) 

      IF (FLOAT(ILG1).LT.FI)                                   THEN
          BOT = (1.-X)*GG(I,J)  +X*GG(1,J) 
          TOP = (1.-X)*GG(I,J+1)+X*GG(1,J+1) 
      ELSE
          BOT = (1.-X)*GG(I,J)  +X*GG(I+1,J) 
          TOP = (1.-X)*GG(I,J+1)+X*GG(I+1,J+1) 
      END IF

      VAL = (1.-Y)*BOT+Y*TOP 

      RETURN
*--------------------------------------------------------------------
 
      END 
      SUBROUTINE gintlm (VAL,F,NI,NJ,FI,FJ,SPVAL)

***    FEB 25/93 - F.MAJAESS (MODIFIED FROM GINTL2 ROUTINE) 

***    LINEAR INTERPOLATION AT POINT (FI,FJ) IN F(NI,NJ) 
***    PROVIDED THE NEEDED GRID POINTS FOR THE COMPUTATION
***    ARE DIFFERENT FROM THE MASK VALUE "SPVAL".

      IMPLICIT none

      INTEGER  NI,NJ, I,J
      REAL     F(NI,NJ),FI,FJ,VAL, SPVAL, X,Y, BOT,TOP

*-------------------------------------------------------------------- 
      VAL = SPVAL
      I   = INT( FI )
      J   = INT( FJ )

      IF (1.0*I.EQ.FI .AND.
     +    1.0*J.EQ.FJ)                                         THEN 
***        GG POINT COINCIDE WITH LL-GRID POINT.
          VAL = F(I,J)
      ELSE
          IF (1.0*I.EQ.FI)                                     THEN 
***            LONGITUDE COINCIDE FOR GG AND LL-GRID.
              IF (J.LT.1 ) J = 1
              IF (J.GE.NJ) J = NJ-1
              Y = FJ-FLOAT( J ) 
              IF (F(I,J  ).NE.SPVAL .AND.
     +            F(I,J+1).NE.SPVAL)
     +            VAL = (1.-Y)*F(I,J) + Y*F(I,J+1)
          ELSE
              IF (1.0*I.EQ.FI)                                 THEN 
***                LATITUDE COINCIDE FOR GG AND LL-GRID.
                  IF (I.LT.1 ) I = 1
                  IF (I.GE.NI) I = NI-1
                  X = FI-FLOAT( I ) 
                  IF (F(I  ,J).NE.SPVAL .AND.
     +                F(I+1,J).NE.SPVAL)
     +                VAL = (1.-X)*F(I,J) + X*F(I+1,J)
              ELSE
***                NEITHER LATITUDES NOR LONGITUDES COINCIDE
                  IF (I.LT.1 ) I = 1
                  IF (I.GE.NI) I = NI-1
                  X = FI-FLOAT( I ) 
                  IF (J.LT.1 ) J = 1 
                  IF (J.GE.NJ) J = NJ-1
                  Y = FJ-FLOAT( J ) 
                  IF (F(I  ,J  ).NE.SPVAL .AND.
     +                F(I+1,J  ).NE.SPVAL .AND.
     +                F(I  ,J+1).NE.SPVAL .AND.
     +                F(I+1,J+1).NE.SPVAL )                    THEN

                      BOT = (1.-X)*F(I,J  )+X*F(I+1,J  )
                      TOP = (1.-X)*F(I,J+1)+X*F(I+1,J+1)

                      VAL = (1.-Y)*BOT+Y*TOP

                  END IF
              END IF
          END IF
      END IF

      RETURN
*----------------------------------------------------------------------- 

      END 
      SUBROUTINE gintc2 (VAL,G,NI,NJ,FI,FJ) 

***   ****   JAN 1975  -  JOHN D. HENDERSON  **** 

***    CUBIC INTERPOLATION AT POINT (FI,FJ) IN GRID G(NI,NJ). 
***    THE GRID MUST BE AT LEAST 4 BY 4 POINTS. 
***    IF (FI,FJ) IS OUTSIDE THE GRID AN EXTRAPOLATION IS DONE. 

      IMPLICIT  none

      INTEGER   NI,NJ, I,J,K,L, IM2,JM2
      REAL      VAL, G(NI,1),AP(4),AQ(4),A(4),
     +          FI,FJ, P,Q, SIXTH, ABSP

*----------------------------------------------------------------------- 
      SIXTH = 1./6. 

***    CALCULATE X-PARAMETERS AND Y-PARAMETERS. 

          I   = INT( FI ) 
      IF (FLOAT(NI).LT.FI)                                     THEN
          I   = NI
      ELSE IF (I.GT.NI-2)                                      THEN
          I   = MAX( 1,NI-2 )
      ELSE IF (I.LT.2)                                         THEN
          I   = MIN( 2,NI )
      END IF
          IM2 = I-2 

      P     =  FI-FLOAT(I) 
      AP(1) = SIXTH*(-P*(P-1.)*(P-2.)) 
      AP(2) =   0.5*(   (P-1.)*(P+1.)*(P-2.)) 
      AP(3) =   0.5*(-P*(P+1.)*(P-2.)) 
      AP(4) = SIXTH*( P*(P+1.)*(P-1.)) 

      ABSP  = ABS( p )

          J   = INT( FJ ) 
      IF (J.GT.NJ-2)                                           THEN
          J   = NJ-2 
      ELSE IF(J.LT.2)                                          THEN
          J   = 2 
      END IF
          JM2 = J-2 

      Q     =  FJ-FLOAT( J ) 
      AQ(1) = SIXTH*(-Q*(Q-1.)*(Q-2.)) 
      AQ(2) =   0.5*(   (Q-1.)*(Q+1.)*(Q-2.)) 
      AQ(3) =   0.5*(-Q*(Q+1.)*(Q-2.)) 
      AQ(4) = SIXTH*( Q*(Q+1.)*(Q-1.)) 

***    INTERPOLATE IN EACH ROW THEN IN THE RESULTING COL FOR VAL. 

      IF (FI.LT.FLOAT(NI) .AND. ABSP.GT.1E-7)                  THEN
          DO  L=1,4 
              A(L) = 0. 
              DO  K=1,4 
                  A(L) = A(L)+AP(K)*G(IM2+K,JM2+L) 
              END DO
          END DO
      ELSE IF (ABSP.LE.1E-7)                                   THEN
          DO  L=1,4
              A(L) = G(I,JM2+L)
          END DO
      ELSE
          DO  L=1,4
              A(L) =      AP(1)*G(IM2+1,JM2+L)
              A(L) = A(L)+AP(2)*G(IM2+2,JM2+L)
              A(L) = A(L)+AP(3)*G(    1,JM2+L)
              A(L) = A(L)+AP(4)*G(    2,JM2+L)
          END DO
      END IF

      VAL = AQ(1)*A(1)+AQ(2)*A(2)+AQ(3)*A(3)+AQ(4)*A(4) 

      RETURN 
*--------------------------------------------------------------------

      END 
      SUBROUTINE gintl2 (VAL,F,NI,NJ,FI,FJ) 

***    LINEAR INTERPOLATION AT POINT (FI,FJ) IN F(NI,NJ). 
***    EXTRAPOLATION IS DONE IF POINT IS OUTSIDE THE GRID 

      IMPLICIT none

      INTEGER  NI,NJ, I,J
      REAL     VAL, F(NI,NJ), FI,FJ, BOT,TOP, X,Y

*-------------------------------------------------------------------- 
          I   = INT( FI ) 
      IF (I.LT.1)                                              THEN
          I   = 1 
      ELSE IF (I.GE.NI)                                        THEN
          I   = MAX( 1,NI-1 )
      END IF

          X   = FI-FLOAT( I ) 

          J   = INT( FJ ) 
      IF (J.LT.1)                                              THEN
          J   = 1  
      ELSE IF (J.GE.NJ)                                        THEN
          J   = NJ-1 
      END IF

          Y   = FJ-FLOAT( J ) 

          IF (NI.GT.1)                                         THEN
              BOT = (1. -X)*F(I,J)   +X*F(I+1,J) 
              TOP = (1. -X)*F(I,J+1) +X*F(I+1,J+1) 
          ELSE
              BOT = F(I,J)
              TOP = F(I,J+1)
          END IF

          VAL = (1.-Y)*BOT+Y*TOP 

      RETURN 
*--------------------------------------------------------------------

      END 
      SUBROUTINE iggsl (GG,ILG1,ILAT,ILEV, GS,LON,ILG,NLAT) 

***    JAN 26/83 - B.DUGAS.

***    INSERTS A 2-D VERTICAL SLICE GS(LON,ILAT) OF WHICH ONLY ILG 
***    POINTS ARE GOOD INTO LATITUDE NLAT OF THE 3-D GAUSSIAN GRID 
***    ARRAY GG(ILG1,ILAT,ILEV). 

***    THE LAST WORD OF EACH ROW IS SET EQUAL TO THE FIRST WORD. 

      IMPLICIT none

      INTEGER  ILG1,ILAT,ILEV,LON,ILG,NLAT,I,L
      REAL     GG(ILG1,ILAT,ILEV),GS(LON,ILEV)
*-------------------------------------------------------------------- 

      DO  L=1,ILEV 
          DO  I=1,ILG
              GG(I,NLAT,L)=GS(I,L)
          END DO
          IF (ILG1.GT.ILG) GG(ILG1,NLAT,L)=GS(1,L) 
      END DO

      RETURN
*--------------------------------------------------------------------

      END 
      SUBROUTINE llfxy (DLAT,DLON,X,Y,D60,DGRW,NHEM)

***   ****   FEB 1975  -  JOHN D. HENDERSON  **** 

***    CALCULATE LATITUDE AND LONGITUDE IN DEGREES OF POINT (X,Y) 
***     MEASURED FROM THE POLE. (LONGITUDE IS POSITIVE EASTWARD). 
***    GRID IS POLAR STEREOGRAPHIC WITH STANDARD LATITUDE AT 60 DEG. 
***     AND GRID SIZE D60 METERS. 

***    ZERO DEGREES LONGITUDE IN THE GRID IS (DGRW) DEGREES 
***     IN MAP COORDINATES. 

***    NHEM 1 = NORTHERN HEMISPHERE.   NHEM 2 = SOUTHERN HEMISPHERE. 
***    1.866025=(1+SIN60),   6.371E+6=EARTH RADIUS IN METERS. 

      IMPLICIT none

      INTEGER  NHEM
      REAL     DLAT,DLON,X,Y,D60,DGRW, RE,RE2,R2, C1

*----------------------------------------------------------------------- 
      RE  = 1.866025*6.371E+6/D60 
      RE2 = RE*RE
      C1  = 180./3.14159 

***    IF POINT IS AT POLE SET COORD TO (0.,90.). 

      DLAT = 90. 
      DLON = 0. 

      IF (X.NE.0. .OR. Y.NE.0.)                                THEN

*         * CALCULATE LONGITUDE IN MAP COORDINATES. 

          IF (X.EQ.0.)                                         THEN
                           DLON = SIGN( 90.,Y ) 
          ELSE
                           DLON = ATAN( Y/X )*C1 
              IF (X.LT.0.) DLON = DLON+SIGN(180.,Y) 
          END IF

*         * ADJUST LONGITUDE FOR GRID ORIENTATION. 

          DLON = DLON-DGRW 
          IF (DLON.GT.+180.)                                   THEN
              DLON = DLON-360. 
          ELSE IF (DLON.LT.-180.)                              THEN
              DLON = DLON+360. 
          END IF

*         * CALCULATE LATITUDE. 

          R2 = X*X+Y*Y 
          DLAT = (RE2-R2)/(RE2+R2) 
          DLAT =  ASIN( DLAT )*C1 

      END IF

***    CHANGE SIGNS IF IN SOUTHERN HEMISPHERE. 

   39 IF (NHEM.EQ.2)                                           THEN
          DLAT = -DLAT 
          DLON=-DLON 
      END IF

      RETURN 
*--------------------------------------------------------------------

      END 
      SUBROUTINE xyfll (X,Y,DLAT,DLON,D60,DGRW,NHEM) 

***   ****   FEB 1975  -  JOHN D. HENDERSON  **** 

***    CALCULATES GRID COORDINATES MEASURED FROM THE POLE OF
***    POINT (DLAT,DLON) GIVEN IN DEGREES. (LONGITUDE IS PO-
***    SITIVE EASTWARD). GRID IS POLAR STEREOGRAPHIC WITH
***    STANDARD LATITUDE AT 60 DEG AND GRID SIZE D60 METERS. 
 
***    ZERO DEGREES LONGITUDE IN THE GRID IS (DGRW) DEGREES 
***     IN MAP COORDINATES. 
***    NHEM 1 = NORTHERN HEMISPHERE.  
***         2 = SOUTHERN HEMISPHERE. 

***    1.866025=(1+SIN60),   6.371E+6=EARTH RADIUS IN METERS. 

      IMPLICIT none

      INTEGER  NHEM
      REAL     DLON,DLAT,D60,DGRW,X,Y
      REAL     GLON,GLAT,RLON,RLAT,SINLAT,R,RE,C2

*----------------------------------------------------------------------- 
      RE = 1.866025*6.371E+6/D60 
      C2 = 3.14159/180. 

                     GLON = DLON 
      IF (NHEM.EQ.2) GLON =-DLON 
                     GLAT = DLAT 
      IF (NHEM.EQ.2) GLAT =-DLAT 

      RLON   = C2*(GLON+DGRW) 
      RLAT   = C2*GLAT 
      SINLAT = SIN( RLAT ) 

      R      = RE*SQRT( (1.-SINLAT)/(1.+SINLAT) ) 
      X      = R*COS( RLON ) 
      Y      = R*SIN( RLON ) 

      RETURN 
*--------------------------------------------------------------------

      END 
      SUBROUTINE lgrdc (DEN,S,NN)

***    NOV 21/77 - J.D.HENDERSON

***    COMPUTES DENOMINATORS FOR LAGRANGIAN CUBIC INTERPOLATOR (LGRIC).
***    NN COORDINATES ARE CONTAINED IN S.
***    DEN(4,NN-3) IS SET TO THE INVERSE OF THE 4 DENOMINATORS
***    FOR EACH OF THE NN-3 POSSIBLE CUBICS.

      IMPLICIT none

      INTEGER  NN,NCUB,NC,NC3,L,K,I
      REAL     DEN(4,1),S(NN),PROD

*--------------------------------------------------------------------
      NCUB=NN-3

      DO 310 NC=1,NCUB
      NC3=NC+3

      L=0
       DO 280 K=NC,NC3
      PROD=1.

      DO 250 I=NC,NC3
      IF(I.NE.K) PROD=PROD*(S(K)-S(I))
  250 CONTINUE

      L=L+1
  280 DEN(L,NC)=1./PROD

  310 CONTINUE
      RETURN
*--------------------------------------------------------------------

      END
      SUBROUTINE lgric (YIN,XC,NX,Y,S,DEN,NN,DL,DR)

***    NOV 29/77 - J.D.HENDERSON

***    PERFORMS LAGRANGIAN INTERPOLATION AT POINTS IN XC(NX)
***     OF THE VARIABLE CONTAINED IN Y(NN), KNOWN AT THE
***    NN COORDINATES CONTAINED IN S(NN).
***    COORDINATES IN S MUST INCREASE TO THE RIGHT (IE 1 TO NN).
***    DEN(4,NN-3) CONTAINS THE INVERTED DENOMINATORS (SEE LGRDC).
***    DL,DR ARE EXTRAPOLATION DERIVATIVES FOR LEFT,RIGHT SIDES.

      IMPLICIT none

      INTEGER  NX,NN,NV,NC,K,I
      REAL     YIN(NX),XC(NX),DL,DR
      REAL     DEN(4,1),S(NN),Y(NN)
      REAL     X,X1,X2,X3,X4,VAL

*--------------------------------------------------------------------
***    LOOP OVER ALL POINTS TO BE INTERPOLATED.

      DO 300 NV=1,NX

          X = XC(NV)

          IF (X.LT.S(1))                                       THEN

***            EXTRAPOLATE TO THE LEFT.

              VAL = Y(1)+DL*(X-S(1))

          ELSE IF (X.GT.S(NN))                                 THEN

***            EXTRAPOLATE TO THE RIGHT.

              VAL = Y(NN)+DR*(X-S(NN))

          ELSE

***            INTERPOLATION:
***            FIRST FIND WHICH CUBIC TO USE.
***            NC IS ALSO THE LEFT POINT OF THE 4 THAT SURROUND X.

              DO 100 I=2,NN
                  K = I
                  IF (X.LT.S(I)) GOTO 200
  100         CONTINUE
  200         NC = K-2

              IF (K.EQ. 2) NC = 1
              IF (K.EQ.NN) NC = NC-1

***            NOW EVALUATE THE LAGRANGIAN POLYNOMIAL AT POINT X.

              X1  = X-S(NC)
              X2  = X-S(NC+1)
              X3  = X-S(NC+2)
              X4  = X-S(NC+3)
              VAL = X2*(X3 *X4)*DEN(1,NC)*Y(NC)
     +            + X1*(X3 *X4)*DEN(2,NC)*Y(NC+1)
     +            +(X1* X2)*X4 *DEN(3,NC)*Y(NC+2)
     +            +(X1* X2)*X3 *DEN(4,NC)*Y(NC+3)

          END IF

          YIN(NV)=VAL

  300 CONTINUE

      RETURN
*--------------------------------------------------------------------

      END
      SUBROUTINE linil (YIN,XC,NX,Y,S,NN,DL,DR)

***    FEB 22/78 - J.D.HENDERSON

***    PERFORMS LINEAR INTERPOLATION AT POINTS IN XC(NX)
***     OF THE VARIABLE CONTAINED IN Y(NN), KNOWN AT THE
***    NN COORDINATES CONTAINED IN S(NN).
***    COORDINATES IN S MUST INCREASE TO THE RIGHT (IE 1 TO NN).
***    DL,DR ARE EXTRAPOLATION DERIVATIVES FOR LEFT,RIGHT SIDES.

      IMPLICIT none

      INTEGER  NX,NN,NV,K,I
      REAL     YIN(NX),XC(NX),X,DL,DR
      REAL     S(NN),Y(NN),VAL,SLOPE

*--------------------------------------------------------------------
***    LOOP OVER ALL POINTS TO BE INTERPOLATED.

      DO 300 NV=1,NX

          X = XC(NV)

          IF (X.LT.S(1))                                          THEN

***            EXTRAPOLATE TO THE LEFT.

              VAL = Y(1)+DL*(X-S(1))

          ELSE IF(X.GT.S(NN))                                     THEN

***            EXTRAPOLATE TO THE RIGHT.

              VAL = Y(NN)+DR*(X-S(NN))

          ELSE

***            INTERPOLATION:
***            FIRST FIND WHICH INTERVAL TO USE.

              DO 100 I=2,NN
                  K = I-1
                  IF (X.LT.S(I)) GO TO 200
  100         CONTINUE

***            NOW PERFORM THE LINEAR INTERPOLATION.

  200         SLOPE = (Y(K+1)-Y(K))/(S(K+1)-S(K))
              VAL   =  Y(K)+SLOPE*(X-S(K))

          END IF

          YIN(NV) = VAL

  300 CONTINUE

      RETURN
*--------------------------------------------------------------------

      END
      SUBROUTINE llcal (DLAT1,DLON1,DLAT2,DLON2,
     +                  DGRW, I1,J1,I2,J2,SLAT, 
     +                  SLON,ILG1,ILAT, NBAD) 

***    MODIFICATION: - NOV 19/04 - B.DUGAS.
***                    NE PLUS CALCULER SLON.
***    MODIFICATION: - AUG 13/90 - M.LAZARE. 
***                    REMOVE HARD-COAT OF SLON BY PASSING IN CALL.

***    AUTHOR:  M.SUTCLIFFE, CCRN - AUG 10/87 

***    GIVEN LAT-LON COORDINATES OF AN AREA, FIND THE INPUT GRID INDICES 
***    CLOSEST TO THESE COORDINATES. ADJUST THE LAT-LON COORDINATES 
***    TO MATCH EXACTLY THE GRID COORDINATES OF THE INDICES. 
***    ALSO CALCULATE DGRW, THE CENTRE LONGITUDE OF THE AREA. 

      IMPLICIT none

      INTEGER  I1,I2,J1,J2,ILG1,ILAT,NBAD
      REAL     DLAT1,DLAT2,DLON1,DLON2,DGRW
      REAL     SLAT(ILAT),SLON(ILG1) 

      INTEGER  J

*-------------------------------------------------------------------------- 

                                             NBAD = 0 
      IF (DLAT1.GE.DLAT2)                    NBAD = 1 
      IF (DLAT1.LT.- 90. .OR. DLAT1.GT. 90.) NBAD = 2 
      IF (DLAT2.LT.- 90. .OR. DLAT2.GT. 90.) NBAD = 3 
      IF (DLON1.LT.-180. .OR. DLON1.GT.180.) NBAD = 4 
      IF (DLON2.LT.-180. .OR. DLON2.GT.180.) NBAD = 5 

      IF (NBAD.NE.0)  RETURN 

      IF (DLON1.LT.0) DLON1 = 360.+DLON1 
      IF (DLON2.LT.0) DLON2 = 360.+DLON2 
                      DLAT1 =  90.+DLAT1 
                      DLAT2 =  90.+DLAT2 

***    FIND THE CORRESPONDING INPUT GRID X CO-ORDINATES 

      I1 = 1 
      I2 = ILG1 
      DO  J=1,ILG1-1 
          IF (SLON(J).LE.DLON1 .AND. SLON(J+1).GE.DLON1) I1 = J
          IF (SLON(J).LE.DLON2 .AND. SLON(J+1).GE.DLON2) I2 = J
      END DO

      IF (SLON(I1+1)-DLON1 .LT. DLON1-SLON(I1)) I1 = I1+1 
      IF (SLON(I2+1)-DLON2 .LT. DLON2-SLON(I2)) I2 = I2+1 

      DLON1 = SLON(I1) 
      DLON2 = SLON(I2) 

***    FIND THE CORRESPONDING INPUT GRID Y INDICES 

      J1 = 1 
      J2 = ILAT 
      DO 300 J=1,ILAT-1 
          IF (SLAT(J).LE.DLAT1 .AND. SLAT(J+1).GT.DLAT1) J1 = J 
          IF (SLAT(J).LE.DLAT2 .AND. SLAT(J+1).GT.DLAT2) J2 = J 
  300 CONTINUE 

      IF (J1.NE. 1  )                                          THEN
      IF (SLAT(J1+1)-DLAT1 .LT. DLAT1-SLAT(J1)) J1 = J1+1 
      END IF
      IF (J2.NE.ILAT)                                          THEN
      IF (SLAT(J2+1)-DLAT2 .LT. DLAT2-SLAT(J2)) J2 = J2+1 
      END IF

      DLAT1 = SLAT(J1) 
      DLAT2 = SLAT(J2) 

***    CALCULATE THE CENTRE LONGITUDE OF THE CHOSEN AREA 

      DGRW = .5*(DLON1+DLON2) 
      IF (I1.GE.I2)                                            THEN 
          IF (DGRW.LT.180.)                                    THEN 
              DGRW = DGRW+180. 
          ELSE 
              DGRW = DGRW-180. 
          END IF 
      END IF

      IF (DGRW.GE.180.) 
     +    DGRW = DGRW-360. 

                         DLAT1 = DLAT1- 90. 
                         DLAT2 = DLAT2- 90. 
      IF (DLON1.GT.180.) DLON1 = DLON1-360. 
      IF (DLON2.GT.180.) DLON2 = DLON2-360. 

      RETURN 
*--------------------------------------------------------------------

      END 
      SUBROUTINE pscal (LY,IP,JP,D60,LX,DGRW,DLAT1,DLON1, 
     +                  DLAT2,DLON2,NHEM,NDEF,NBAD) 

***    CCRN AUG 10/87. M.SUTCLIFFE, R.DALEY. 

***    CALCULATES PARAMETERS LY,IP,JP,D60, FOR POLAR STEREOGRAPHIC 
***    CODES USING NUMBER OF X GRIDPOINTS (LX), GREENWICH ANGLE (DGRW), 
***    AND LAT-LONG OF LOWER LEFT HAND CORNER (DLAT1,DLON1) AND 
***    UPPER RIGHT HAND CORNER (DLAT2,DLON2).  (BECAUSE LX,LY,IP,JP 
***    MUST ALL BE INTEGERS, THE ACTUAL AREA COVERED WILL BE SLIGHTLY 
***    LARGER THAN SPECIFIED).  NEW VALUES OF DLAT1,DLON1,DLAT2,DLON2 
***    ARE CALCULATED AS WELL. 

***    NHEM REPRESENTS THE (N,S) POLAR STEREOGRAPHIC PROJECTION TO BE USED 
***    GIVEN VALUES OF (1,2). 

***    IF THE DEFAULT (NDEF=0) IS GIVEN, THE OLD CCRN STANDARD PROJECTION 
***    AREA IS RE-CALCULATED EXACTLY, BY REDEFINING THE VALUES OF E AND IOFF. 

      IMPLICIT none

      INTEGER  LY,IP,JP,LX,NHEM,NDEF,NBAD
      REAL     DGRW,DLAT1,DLON1,DLAT2,DLON2,D60

      INTEGER  IOFF
      REAL     E,X1,X2,Y1,Y2,DIST

      EXTERNAL XYFLL,LLFXY

      DATA     E    / 1.0 / 
      DATA     IOFF / 1   / 

*------------------------------------------------------------------------ 
                                             NBAD = 0 
      IF (DLAT1.LT.- 90. .OR. DLAT1.GT. 90.) NBAD = 2 
      IF (DLAT2.LT.- 90. .OR. DLAT2.GT. 90.) NBAD = 3 
      IF (DLON1.LT.-180. .OR. DLON1.GT.180.) NBAD = 4 
      IF (DLON2.LT.-180. .OR. DLON2.GT.180.) NBAD = 5 

      IF (NBAD.NE.0) RETURN 

      IF (NDEF.EQ.0)                                           THEN 
          E    = 1.E-2 
          IOFF = 0 
      END IF 

***    FIND GRID LENGTH (D60). 

      CALL XYFLL( X1,Y1, DLAT1,DLON1,1.E5,DGRW, NHEM ) 
      CALL XYFLL( X2,Y2, DLAT2,DLON2,1.E5,DGRW, NHEM ) 
      DIST = 1.E5*ABS( X1-X2 ) 
      D60 = DIST/FLOAT( LX-1 - IOFF ) 

***    CALCULATE POLE POSITIONS (IP,JP) 

      CALL XYFLL( X1,Y1, DLAT1,DLON1,D60,DGRW, NHEM ) 
      IP = -X1 + 1. + E 
      JP = -Y1 + 1. + E 
      X1 = -FLOAT( IP-1 ) 
      Y1 = -FLOAT( JP-1 ) 
      CALL LLFXY( DLAT1,DLON1, X1,Y1,D60,DGRW, NHEM ) 

***    CALCULATE NUMBER OF GRIDPOINTS IN Y DIRECTION (LY) 

      CALL XYFLL( X2,Y2, DLAT2,DLON2,D60,DGRW, NHEM ) 
      CALL XYFLL( X1,Y1, DLAT1,DLON1,D60,DGRW, NHEM ) 
      LY =  ABS( Y2-Y1 ) + 1. + E 
      X2 = X1 + FLOAT( LX-1 ) 
      Y2 = Y1 + FLOAT( LY-1 ) 
      CALL LLFXY( DLAT2,DLON2, X2,Y2,D60,DGRW, NHEM ) 

      RETURN 
*-------------------------------------------------------------------- 

      END 
      SUBROUTINE llegg2 (GG,ILG1,ILAT, DLAT, GLL,NLG1,NLAT) 

***    DEC  3/79 - J.D.HENDERSON 

***    EXTRACTS A GLOBAL GAUSSIAN GRID GG(ILG1,ILAT) FROM A 
***    GLOBAL LAT-LONG GRID GLL(NLG1,NLAT). EACH GG POINT
***    TAKES THE VALUE OF THE CLOSEST POINT IN GLL. 

***    GG HAS ILG EQUALLY SPACED LONGITUDES. ILAT IS EVEN AND THE 
***    EQUATOR AND POLES ARE NOT ONE OF THE LATITUDE CIRCLES. 
***    DLAT(ILAT) CONTAINS GG LATITUDES IN DEGREES S. TO N. 

***    THE LAT-LONG GRID GLL HAS NLG EQUALLY SPACED LONGITUDES 
***    AND NLAT EQUALLY SPACED LATITUDES. THE POLES ARE INCLUDED
***    AS THE TOP AND BOTTOM ROWS. 

***    FOR BOTH GRIDS THE LEFT SIDE ARE REPEATED ON THE RIGHT
***    SIDE WHEN ILG1 OR NLG1 ARE UN-EVEN NUMBERS.

      IMPLICIT none

      INTEGER     ILG1,ILAT,     NLG1,NLAT
      REAL     GG(ILG1,ILAT),GLL(NLG1,NLAT) 
      REAL     DLAT(ILAT) 

      REAL*8   DLON,IDEG
      INTEGER  ILG,NLG,   I,J,  II,JJ
      REAL     DEGX,DEGY, DROW, DX,DY, X,Y

*----------------------------------------------------------------------- 
                              ILG = ILG1-1
      IF (MOD( ILG1,2 ).EQ.0) ILG = ILG1
                              NLG = NLG1-1
      IF (MOD( NLG1,2 ).EQ.0) NLG = NLG1

      DEGX = 360./FLOAT( NLG ) 
      DEGY = 180./FLOAT( NLAT-1 ) 
      IDEG = 360./DBLE( ILG )

***    THE ROWS ARE DONE FROM SOUTH TO NORTH. 
***    DROW = DEGREES OF ROW J FROM THE SOUTH POLE. 

      DO 240 J=1,ILAT 

          DROW = DLAT(J)+90. 
          Y    = DROW/DEGY+1. 
          JJ   = INT( Y ) 
          DY   = Y-FLOAT( JJ ) 
          IF (DY.GT.0.5) 
     +    JJ   = JJ+1 

***    EXTRACT THE POINTS FROM LEFT TO RIGHT IN EACH ROW. 
***    DLON = DEGREES OF POINT I EASTWARD FROM GREENWICH. 

          DLON = 0.0

          DO 220 I=1,ILG 

              X    = DLON/DEGX+1. 
              DLON = DLON+IDEG
              II   = INT( X ) 
              DX   = X-FLOAT( II ) 

              IF (DX.GT.0.5)
     +        II   = II+1

              IF (II.GT.NLG)
     +        II   = 1

              GG(I,J)=GLL(II,JJ) 

  220     CONTINUE 

          IF (ILG.NE.ILG1) GG(ILG1,J) = GG(1,J)

  240 CONTINUE 

      RETURN 
*-------------------------------------------------------------------- 

      END 
      SUBROUTINE lliggm (GG,ILG1,ILAT,DLAT, GLL,NLG1,NLAT,
     +                                    ALON,ALAT,SPVAL)

***    FEB 25/93 - F.MAJAESS (MODIFIED FROM LLIGG2 ROUTINE) 

***    LINEARILY INTERPOLATES GLOBAL GAUSSIAN GRID 
***    GG(ILG1,ILAT) FROM GLOBAL LAT-LONG GRID 
***    GLL(NLG1,NLAT) AT THE LONGITUDES IN ALON(NLG1) AND 
***    LATITUDES IN ALAT(NLAT).
***    VALUES EQUAL TO "SPVAL" ARE TREATED AS BAD/MISSING
***    DATA.
***    DLAT = LAT  (DEG) OF GG  ROWS.
***    ALON = LONG (DEG) OF GLL COLUMNS.
***    ALAT = LAT  (DEG) OF GLL ROWS.

      IMPLICIT none

      INTEGER  ILG1,ILAT,NLG1,NLAT
      REAL     GG(ILG1,ILAT),GLL(NLG1,NLAT) 
      REAL     DLAT(ILAT),ALON(NLG1),ALAT(NLAT),SPVAL

      REAL     X,Y,DROW,DLON,VAL
      INTEGER  ILG,NLG,I,J,II,JJ

      EXTERNAL GINTLM

*-------------------------------------------------------------------- 
                              ILG = ILG1-1
      IF (MOD( ILG1,2 ).EQ.0
     +.OR.     ILG1    .LE.1) ILG = ILG1 
                              NLG = NLG1-1
      IF (MOD( NLG1,2 ).EQ.0
     +.OR.     NLG1    .LE.1) NLG = NLG1 

***    ROWS ARE DONE FROM SOUTH TO NORTH.
***    DROW = DEGREES OF ROW J FROM S POLE.
***    (X,Y) = COORD OF GG(I,J) IN LAT-LONG GRID.

      DO  400 J=1,ILAT 

          DROW=DLAT(J)
          IF (DROW.GE.ALAT( 1  ) .AND.
     +        DROW.LE.ALAT(NLAT))                              THEN
              Y = 0.0
              DO  JJ=1,NLAT-1
                  IF (DROW.GE.ALAT(JJ  ) .AND.
     +                DROW.LE.ALAT(JJ+1))                      THEN
                      Y = JJ
     +                  + (DROW-ALAT(JJ))/(ALAT(JJ+1)-ALAT(JJ))
                      GOTO 100
                  END IF
              END DO
          ELSE
	      IF (DROW.LT.ALAT(1))                             THEN
                  Y = 1.0
              ELSE
                  Y = FLOAT( NLAT )
              END IF
          END IF

  100     CONTINUE

          IF (NLG.GT.1 .AND. ILG.GT.1)                         THEN

***            POINTS ARE INTERPOLATED FROM LEFT TO RIGHT. 
***            DLON = DEGREES OF POINT I EAST FROM GREENWICH.

              DO  300 I=1,ILG

                  DLON = FLOAT( I-1 )/FLOAT( ILG )*360. 

                  IF (DLON.GE.ALON( 1 ) .AND.
     +                DLON.LE.ALON(NLG))                       THEN
                      X = 0.0
                      DO  II=1,NLG
                          IF (DLON.GE.ALON(II  ) .AND.
     +                        DLON.LE.ALON(II+1))              THEN
                              X = II
     +                          +  ( DLON - ALON(II))
     +                          / (ALON(II+1)-ALON(II))
                              GOTO 200
                          END IF
                      END DO
                  ELSE
                      IF (DLON.LT.ALON(NLG)) DLON=DLON+360. 
                      X = NLG+(DLON-ALON(NLG))/(360.-ALON(NLG))
                  END IF

 200              CONTINUE

                  VAL     = SPVAL
                  CALL GINTLM( VAL,GLL,NLG1,NLAT,X,Y, SPVAL )
                  GG(I,J) = VAL 

 300          CONTINUE

***            REPEAT GREENWHICH MERIDIAN AS EXTRA CYCLIC LONGITUDE. 

              IF (MOD(ILG1,2).NE.0) GG(ILG1,J) = GG(1,J)

          ELSE

***            WE HAVE A ZONAL ARRAY.

              VAL     = SPVAL
              CALL GINTLM( VAL,GLL,1,NLAT,1.0,Y, SPVAL )
              GG(1,J) = VAL 

          END IF

  400 CONTINUE

      RETURN
*-------------------------------------------------------------------- 

      END 
      SUBROUTINE lligg2 (GG,ILG1,ILAT, DLAT, GLL,NLG1,NLAT, INTERP)

***    APR 30/80 - J.D.HENDERSON 

***    INTERPOLATES GLOBAL GAUSSIAN GRID GG(ILG1,ILAT) 
***    FROM GLOBAL LAT-LONG GRID GLL(NLG1,NLAT). 
***    DLAT = LAT (DEG) OF GG ROWS. 
***    INTERP = (1,3) FOR (LINEAR,CUBIC) INTERPOLATION. 
***             (OTHERWISE THE GRID GG IS SET TO ZERO). 

      IMPLICIT none

      INTEGER  ILG1,ILAT, NLG1,NLAT, INTERP
      REAL     GG(ILG1,ILAT),GLL(NLG1,NLAT) 
      REAL     DLAT(ILAT)

      EXTERNAL GINTL2,GINTC2

      REAL*8   DLON,IDEG
      INTEGER  ILG,NLG, I,J
      REAL     DEGX,DEGY, VAL, DROW, X,Y

*-------------------------------------------------------------------- 
                              ILG = ILG1-1
      IF (MOD( ILG1,2 ).EQ.0
     +.OR.     ILG1    .LE.1) ILG = ILG1 
                              NLG = NLG1-1
      IF (MOD( NLG1,2 ).EQ.0
     +.OR.     NLG1    .LE.1) NLG = NLG1 

      DEGX = 360./FLOAT( NLG ) 
      DEGY = 180./FLOAT( NLAT-1 ) 
      IDEG = 360./DBLE( ILG )

***    ROWS ARE DONE FROM SOUTH TO NORTH. 
***    DROW = DEGREES OF ROW J FROM S POLE. 
***    (X,Y) = COORD OF GG(I,J) IN LAT-LONG GRID. 

      DO 240 J=1,ILAT 

          DROW = DLAT(J)+90. 
          Y    = DROW/DEGY+1. 

***    POINTS ARE INTERPOLATED FROM LEFT TO RIGHT. 
***    DLON = DEGREES OF POINT I EAST FROM GREENWICH. 

          DLON = 0.0

          DO 220 I=1,ILG 

              X    = DLON/DEGX+1. 
              DLON = DLON+IDEG

              IF (INT( X-NLG ).GT.0.5)
     +        X    = FLOAT( NLG+1 )-X

              VAL=0. 

              IF (INTERP.EQ.1) CALL GINTL2( VAL,GLL,NLG1,NLAT,X,Y ) 
              IF (INTERP.EQ.3) CALL GINTC2( VAL,GLL,NLG1,NLAT,X,Y ) 

              GG(I,J) = VAL 

  220     CONTINUE 

          IF (ILG.NE.ILG1) GG(ILG1,J) = GG(1,J)

  240 CONTINUE 

      RETURN 
*-------------------------------------------------------------------- 

      END 
      SUBROUTINE gcround (GCROW, IS,IF) 

***    JUL 25/84 - R.LAPRISE. 

***    ROUND VALUES OF GCROW TO NEAREST ELEMENT OF SET (-1., 0., 1.). 

      IMPLICIT none

      INTEGER  IS,IF,I
      REAL     GCROW(IF),GC

*---------------------------------------------------------------------- 
      DO 100 I=IS,IF 

                               GC = 0.0
          IF (GCROW(I).LT.-.5) GC =-1.0 
          IF (GCROW(I).GT.+.5) GC =+1.0 

          GCROW(I) = GC 

  100 CONTINUE 

      RETURN 
*-------------------------------------------------------------------- 

      END 
      SUBROUTINE wheneq (N,ARRAY,INC,TARGET,INDEX,NVAL)

      IMPLICIT none

      INTEGER  N,INDEX(N),INC,NVAL, I,INA
      REAL      ARRAY(N), TARGET
      INTEGER  IARRAY(N),ITARGET

*---------------------------------------------------------------------- 
      NVAL=0

                    INA = 1
      IF (INC.LT.0) INA = (-INC)*(N-1)+1

      DO  I=1,N
          IF (ARRAY(INA).EQ.TARGET)                            THEN
              NVAL         = NVAL+1
              INDEX(NVAL ) = I
          END IF
          INA = INA+INC
      END DO

      RETURN
*-------------------------------------------------------------------- 

      ENTRY wheneqi (N,IARRAY,INC,ITARGET,INDEX,NVAL)

      NVAL=0

                    INA = 1
      IF (INC.LT.0) INA = (-INC)*(N-1)+1

      DO  I=1,N
          IF (IARRAY(INA).EQ.ITARGET)                            THEN
              NVAL         = NVAL+1
              INDEX(NVAL ) = I
          END IF
          INA = INA+INC
      END DO

      RETURN

*-------------------------------------------------------------------- 
      END
      SUBROUTINE hralr (GLR,ILG1,ILAT,DLON,DLAT, 
     +                  GHR,NLG1,NLAT,GCLR,GCHR,ICHOICE,OK,
     +                  LONBAD,LATBAD)

***    JUNE 10/86 - M.LAZARE 

***    TRANSFORMS DATA ON A HIGH RESOLUTION GRID(NLG1 X NLAT) TO A LOWER 
***    RESOLUTION GRID(ILG1 X ILAT), BY TAKING INTO ACCOUNT ALL THE HIGH- 
***    RESOLUTION GRID POINTS IN EACH LOWER-RESOLUTION GRID SQUARE HAVING THE 
***    SAME GROUND COVER AS THAT OF THE LOW-RESOLUTION GRID POINT. 

***    THE LONGITUDES AND LATITUDES OF THE HIGH-RESOLUTION FIELD GHR ARE 
***    CONTAINED IN THE ARRAYS RLON AND RLAT RESPECTIVELY, WHILE THOSE OF THE 
***    LOW-RESOLUTION FIELD GLR ARE CONTAINED IN DLON AND DLAT RESPECTIVELY. 

***    GCHR IS THE HIGH-RESOLUTION GROUND COVER ARRAY(NLG1XNLAT) AND GCLR IS 
***    THE LOW-RESOLUTION GROUND COVER ARRAY(ILG1XILAT). 

***    ONLY USE THIS ROUTINE FOR SURFACE FIELDS WHOSE RESOLUTION IS THE 
***    SAME AS THE HIGH-RESOLUTION GROUND COVER FIELD. OTHERWISE, USE THE 
***    LESS RESTRICTIVE VERSION CALLED HRTOLR WHICH DOES NOT CONSIDER THE 
***    GROUND COVER IN PERFORMING THE TRANSFORMATION. 

***    IF ICHOICE.EQ.0: THE MOST FREQUENT VALUE IN THE SQUARE IS RETURNED. IF 
***                 THERE IS A TIE, THE VALUE OF THE POINT CLOSEST TO THE 
***                 LOW-RESOLUTION GRID POINT IS RETURNED. 
***    IF ICHOICE.EQ.1: THE PARTIAL BOX AREA-AVERAGED VALUE IS RETURNED. 

***    IF NO HIGH-RESOLUTION POINTS ARE FOUND WITHIN THE LOW-RESOLUTION GRID 
***    SQUARE, THE SUBROUTINE RETURNS WITH OK=.FALSE. AND PASSES BACK THE 
***    LOCATION OF THE BAD POINT IN (I,J) COORDINATES DEFINED BY 
***    (LONBAD,LATBAD). 

      IMPLICIT    none

      INTEGER     MAXNLG
      PARAMETER ( MAXNLG = 1999 )
      INTEGER     MAXNLA
      PARAMETER ( MAXNLA = 999  )

      LOGICAL     OK 
      INTEGER     ILG1,ILAT,NLG1,NLAT,LONBAD,LATBAD,ICHOICE
      REAL        GLR(ILG1,ILAT),GCLR(ILG1,ILAT),DLAT(ILAT)
      REAL        GHR(NLG1,NLAT),GCHR(NLG1,NLAT),DLON(ILG1) 

      REAL        LATM,LATP,LONM,LONP,DEGX,DEGY,TEMP,PI
      REAL        FACTLAT,FACTL,SUM,SLATP,SLATM,SLONM,SLONP
      REAL        AREA,BIGAREA,DISTRY,DISTLON,DISTMIN,DISTLAT
      REAL        RLAT(MAXNLA),RLON(MAXNLG)
      REAL        VAL(MAXNLG), DIST(MAXNLG)

      INTEGER     I,J,II,JJ,LL,IJK,N,NLG,ILG,NLPTS,NPTS
      INTEGER     NDESERT,NEWMAX,NT,NUM,NVAL,NPP,NCHOICE
      INTEGER     LOCAT(MAXNLG), LATL(MAXNLA),LATH(MAXNLA)
      INTEGER     LOCFST(MAXNLG),IVAL(MAXNLG),NMAX(MAXNLG)
      INTEGER     LONL(MAXNLG),  LONH(MAXNLG)

*------------------------------------------------------------------------------ 
      PI = 4.0*ATAN( 1.D0 )

      LONBAD=0 
      LATBAD=0 
      OK=.TRUE. 

                              ILG = ILG1-1
      IF (MOD( ILG1,2 ).EQ.0) ILG = ILG1 
                              NLG = NLG1-1
      IF (MOD( NLG1,2 ).EQ.0) NLG = NLG1 

***    DEFINE LONGITUDE AND LATITUDE VECTORS FOR HIGH-RESOLUTION GRID. 

      DEGX=360./(FLOAT(NLG1-1)) 
      DEGY=180./(FLOAT(NLAT-1)) 

      DO  10 I=1,NLG1 
   10 RLON(I)=DEGX*FLOAT(I-1) 
      DO  20 J=1,NLAT 
   20 RLAT(J)=-90.+DEGY*FLOAT(J-1) 

***    FIND HIGH-RESOLUTION GRID POINTS WITHIN EACH LOW-RESOLUTION GRID SQUARE. 
***    THIS REQUIRES SPECIAL HANDLING OF POINTS NEAR THE POLES OR AT THE 
***    GREENWICH MERIDIAN. 

      DO 40 J=1,ILAT 
        IF(J.EQ.ILAT) THEN 
          LATP=90. 
        ELSE 
          LATP=0.5*(DLAT(J+1)+DLAT(J)) 
        ENDIF 
        IF(J.EQ.1) THEN 
          LATM=-90. 
        ELSE 
          LATM=0.5*(DLAT(J-1)+DLAT(J)) 
        ENDIF 
        DO 30 JJ=1,NLAT 
          IF(JJ.NE.1) THEN 
            IF(RLAT(JJ-1).LT.LATM.AND.RLAT(JJ).GE.LATM) LATL(J)=JJ 
          ELSE 
            IF(RLAT(JJ).GE.LATM) LATL(J)=JJ 
          ENDIF 
          IF(JJ.NE.NLAT) THEN 
            IF(RLAT(JJ).LE.LATP.AND.RLAT(JJ+1).GT.LATP) LATH(J)=JJ 
          ELSE 
            IF(RLAT(JJ).LE.LATP) LATH(J)=JJ 
          ENDIF 
   30   CONTINUE 
   40 CONTINUE 

      DO 75 I=1,ILG 
        IF(I.EQ.1) THEN 
          LONM=0.5*(DLON(ILG)+360.) 
          LONP=0.5*(DLON(I+1)+DLON(I)) 
          IF(LONM.GT.RLON(NLG)) LONM=0. 
        ELSE IF (I.EQ.ILG) THEN
          LONM=0.5*(DLON(I-1)+DLON(I)) 
          LONP=0.5*(DLON(ILG)+360.)      
        ELSE
          LONM=0.5*(DLON(I-1)+DLON(I)) 
          LONP=0.5*(DLON(I+1)+DLON(I)) 
        ENDIF 
        IF(RLON(1).GE.LONM)                           LONL(I)=1
        IF(RLON(1).LE.LONP.AND.RLON(2).GT.LONP)       LONH(I)=1
        DO 50 II=2,NLG-1
          IF(RLON(II).GE.LONM.AND.RLON(II-1).LT.LONM) LONL(I)=II 
          IF(RLON(II).LE.LONP.AND.RLON(II+1).GT.LONP) LONH(I)=II 
   50   CONTINUE 
        IF(RLON(NLG).GE.LONM.AND.RLON(NLG-1).LT.LONM) LONL(I)=NLG
        IF(RLON(NLG).LE.LONP)                         LONH(I)=NLG
        IF(LONL(I).GT.LONH(I)) LONH(I)=LONH(I)+NLG 
   75 CONTINUE 

***    LOOP OVER POINTS IN THE LOW-RESOLUTION FIELD, KEEPING TRACK OF ALL 
***    HIGH-RESOLUTION POINTS WITHIN EACH LOW-RESOLUTION GRID SQUARE 
***    HAVING THE SAME GROUND COVER AS THE LOW-RESOLUTION GRID SQUARE. 

      DO 900 J=1,ILAT 
        FACTL=(COS(DLAT(J)*PI/180.))**2 
        DO 800 I=1,ILG 
          NDESERT=0 
          NLPTS=0 
          IF(ICHOICE.EQ.1) THEN 

***          FOR CONTINUOUS FIELDS, AREA-AVERAGE THESE POINTS. 

            BIGAREA=0. 
            SUM=0. 
            DO 200 JJ=LATL(J),LATH(J) 
              IF(JJ.EQ.NLAT) THEN 
                SLATP=90. 
              ELSE 
                SLATP=0.5*(RLAT(JJ+1)+RLAT(JJ)) 
              ENDIF 
              IF(JJ.EQ.1) THEN 
                SLATM=-90. 
              ELSE 
                SLATM=0.5*(RLAT(JJ-1)+RLAT(JJ)) 
              ENDIF 
              FACTLAT=(COS(RLAT(JJ)*PI/180.))*(SLATP-SLATM)*PI/180. 
              DO 100 LL=LONL(I),LONH(I) 
                IF(LL.GT.NLG) THEN 
                  II=LL-NLG 
                ELSE 
                  II=LL 
                ENDIF 
                IF(GCHR(II,JJ).EQ.999.) NDESERT=NDESERT+1 
                IF(GCHR(II,JJ).EQ.GCLR(I,J)) THEN 
                  NLPTS=NLPTS+1 
                  IF(II.EQ.1) THEN 
                    SLONM=0.5*(RLON(NLG)+360.)-360. 
                    SLONP=0.5*(RLON(II+1)+RLON(II)) 
                  ELSE IF(II.EQ.NLG) THEN
                    SLONM=0.5*(RLON(II-1)+RLON(II)) 
                    SLONP=0.5*(RLON(NLG)+360.)
                  ELSE 
                    SLONM=0.5*(RLON(II-1)+RLON(II)) 
                    SLONP=0.5*(RLON(II+1)+RLON(II)) 
                  ENDIF 
                  AREA=FACTLAT*((SLONP-SLONM)*PI/180.) 
                  SUM=SUM+AREA*GHR(II,JJ) 
                  BIGAREA=BIGAREA+AREA 
                ENDIF 
  100         CONTINUE 
  200       CONTINUE 
            IF(NLPTS.NE.0) THEN 
              GLR(I,J)=SUM/BIGAREA 
            ENDIF 
          ELSE 

***          FOR DISCRETE-VALUED FIELDS, KEEP TRACK OF EACH QUALIFYING HIGH- 
***          RESOLUTION POINT VALUE IN ARRAY VAL, AND THEIR DISTANCE FROM THE 
***          LOW-RESOLUTION GRID POINT IN ARRAY DIST. 

            DO 400 JJ=LATL(J),LATH(J) 
              DISTLAT=DLAT(J)-RLAT(JJ) 
              DO 300 LL=LONL(I),LONH(I) 
                IF(LL.GT.NLG) THEN 
                  II=LL-NLG 
                ELSE 
                  II=LL 
                ENDIF 
                IF(GCHR(II,JJ).EQ.GCLR(I,J)) THEN 
                  NLPTS=NLPTS+1 
                  VAL(NLPTS)=GHR(II,JJ) 
                  DISTRY=ABS(DLON(I)-RLON(II)) 
                  DISTLON=AMIN1(DISTRY,(360.-DISTRY)) 
                  DIST(NLPTS)=(FACTL*(DISTLON**2)+DISTLAT**2)* 
     +                        ((PI/180.)**2) 
                ENDIF 
  300         CONTINUE 
  400       CONTINUE 
          ENDIF 

***        RETURN WITH OK=.FALSE. IF NO HIGH-RESOLUTION POINTS IN 
***        LOW-RESOLUTION GRID SQUARE DEFINED BY (LONBAD,LATBAD)=(I,J). 
***        FOR SPECIAL CASE OF CANOPY ALBEDOES WHERE HIGH-RESOLUTION 
***        GROUND COVER FIELD ARE ARTIFICIALLY SET TO 999 IN MAIN PROGRAM 
***        IF HIGH-RESOLUTION CANOPY ALBEDO IS ZERO (I.E. PURE DESERT), 
***        ASSIGN VALUE OF ZERO TO RESULTING LOW-RESOLUTION ALBEDO FIELD 
***        IF ALL HIGH-RESOLUTION GRID POINTS ARE PURE DESERT. THIS NEEDS 
***        TO BE DONE EXPLICITLY HERE BECAUSE LOW-RESOLUTION GROUND COVER 
***        FIELD CANNOT BE ASSIGNED VALUES OF 999 AND HENCE WILL GET 
***        NLPTS=0 IF LOW-RESOLUTION GRID SQUARE IS COMPLETELY DESERT. 

          NPTS=(LATH(J)-LATL(J)+1)*(LONH(I)-LONL(I)+1) 
          IF(NLPTS.EQ.0) THEN 
            IF(NDESERT.EQ.NPTS.AND.GCLR(I,J).EQ.+1.) THEN 
               GLR(I,J)=0. 
               NLPTS=NPTS 
            ELSE 
               OK=.FALSE. 
               LONBAD=I 
               LATBAD=J 
               RETURN 
            ENDIF 
          ENDIF 

***        IF ONLY ONE HIGH-RESOLUTION GRID POINT EXISTS, RETURN ITS VALUE. 

          IF(ICHOICE.EQ.0.AND.NLPTS.EQ.1) THEN 
            GLR(I,J)=VAL(1) 
          ELSE IF(ICHOICE.EQ.0.AND.NLPTS.GT.1) THEN 

***          CALCULATE NUMBER OF DISTINCT-VALUED POINTS (DEFINED BY NT), THEIR 
***          FIRST-OCCURRING LOCATION (DEFINED BY ARRAY LOCFST) AND THE 
***          ASSOCIATED LOCATIONS WHERE THEY OCCUR IN ARRAY VAL (DEFINED BY 
***          ARRAY NMAX). 

            NEWMAX=0 
            DO 500 N=1,NLPTS 
              IF(N.EQ.1) THEN 
                NT=0 
                NUM=0 
              ELSE 
                CALL WHENEQ(N-1,VAL,1,VAL(N),IVAL,NUM) 
              ENDIF 
              IF(NUM.EQ.0) THEN 
                NT=NT+1 
                CALL WHENEQ(NLPTS,VAL,1,VAL(N),IVAL,NVAL) 
                LOCFST(NT)=N 
                NMAX(NT)=NVAL 
                IF(NVAL.GT.NEWMAX) NEWMAX=NVAL 
              ENDIF 
  500       CONTINUE 

***          IF ONLY ONE VALUE EXISTS IN THE LOW-RESOLUTION GRID SQUARE, 
***          RETURN ITS VALUE. 

            IF(NT.EQ.1) THEN 
              GLR(I,J)=VAL(LOCFST(1)) 
            ELSE 

***            IF ONLY ONE UNIQUE-VALUED POINT HAS THE MOST FREQUENT OCCURRENCE 
***            NEWMAX, RETURN ITS VALUE. 

              CALL WHENEQI(NT,NMAX,1,NEWMAX,IVAL,NVAL) 
              IF(NVAL.EQ.1) THEN 
                GLR(I,J)=VAL(LOCFST(IVAL(1))) 
              ELSE 

***              OTHERWISE, SEARCH THROUGH THE UNIQUE-VALUED POINTS TO 
***              DETERMINE THE VALUES HAVING THE MOST FREQUENT OCCURRENCE AS 
***              NEWMAX AND THEIR LOCATION IN ARRAY VAL (DEFINED BY ARRAY 
***              LOCAT). 

                NPP=0 
                DO 600 N=1,NT 
                  IF(NMAX(N).EQ.NEWMAX) THEN 
                    TEMP=VAL(LOCFST(N)) 
                    CALL WHENEQ(NLPTS,VAL,1,TEMP,IVAL,NVAL) 
                    DO 550 IJK=1,NVAL 
                      NPP=NPP+1 
                      LOCAT(NPP)=IVAL(IJK) 
  550               CONTINUE 
                  ENDIF 
  600           CONTINUE 

***              CHOOSE THE CLOSEST POINT TO THE TARGET HAVING THAT MAXIMUM- 
***              OCCURRING VALUE. 

                DISTMIN=1.E20 
                DO 700 N=1,NPP 
                  IF(DIST(LOCAT(N)).LT.DISTMIN) THEN 
                    DISTMIN=DIST(LOCAT(N)) 
                    NCHOICE=LOCAT(N) 
                  ENDIF 
  700           CONTINUE 
                GLR(I,J)=VAL(NCHOICE) 
              ENDIF 
            ENDIF 
          ENDIF 
  800   CONTINUE 

***      REPEAT GREENWHICH MERIDIAN AS EXTRA CYCLIC LONGITUDE. 

        IF (MOD(ILG1,2).NE.0) GLR(ILG1,J) = GLR(1,J)

  900 CONTINUE 

      RETURN 
*------------------------------------------------------------------------------ 

      END 
      SUBROUTINE hrtolr (GLR,ILG1,ILAT,DLON,DLAT, 
     +                   GHR,NLG1,NLAT,ICHOICE,OK,
     +                   LONBAD,LATBAD) 

***    JUNE 10/86 - M.LAZARE 

***    TRANSFORMS DATA ON A HIGH RESOLUTION GRID(NLG1 X NLAT) TO A LOWER 
***    RESOLUTION GRID(ILG1 X ILAT), BY TAKING INTO ACCOUNT ALL THE HIGH- 
***    RESOLUTION GRID POINTS IN EACH LOWER-RESOLUTION GRID SQUARE. 

***    THE LONGITUDES AND LATITUDES OF THE HIGH-RESOLUTION FIELD GHR ARE 
***    CONTAINED IN THE ARRAYS RLON AND RLAT RESPECTIVELY, WHILE THOSE OF THE 
***    LOW-RESOLUTION FIELD GLR ARE CONTAINED IN DLON AND DLAT RESPECTIVELY. 

***    DON'T USE THIS ROUTINE FOR SURFACE FIELDS WHOSE RESOLUTION IS THE 
***    SAME AS THE HIGH-RESOLUTION GROUND COVER FIELD. INSTEAD, USE THE 
***    SUBROUTINE HRALR. 

***    IF ICHOICE.EQ.0: THE MOST FREQUENT VALUE IN THE SQUARE IS RETURNED. IF 
***                 THERE IS A TIE, THE VALUE OF THE POINT CLOSEST TO THE 
***                 LOW-RESOLUTION GRID POINT IS RETURNED. 
***    IF ICHOICE.EQ.1: THE PARTIAL BOX AREA-AVERAGED VALUE IS RETURNED. 
***    IF ICHOICE.EQ.2: THE FIELD BEING SMOOTHED IS THE GROUND COVER FIELD. 
***                 THIS IS HANDLED AS IN ICHOICE=0, EXCEPT THAT A DECISION 
***                 IS MADE FIRST ON THE MOST FREQUENT VALUE BEING LAND OR 
***                 WATER (FROZEN OR OPEN). FOR RESULTING NON-LAND POINTS, 
***                 A FINAL DECISION IS MADE ON WHETHER THE MOST FREQUENT 
***                 VALUE IS WATER OR SEA-ICE. 

***    IF NO HIGH-RESOLUTION POINTS ARE FOUND WITHIN THE LOW-RESOLUTION GRID 
***    SQUARE, THE SUBROUTINE RETURNS WITH OK=.FALSE. AND PASSES BACK THE 
***    LOCATION OF THE BAD POINT IN (I,J) COORDINATES DEFINED BY 
***    (LONBAD,LATBAD). 

      IMPLICIT    REAL (A-H,O-Z), INTEGER (I-N)

      INTEGER     MAXNLG
      PARAMETER ( MAXNLG = 1999 )
      INTEGER     MAXNLA
      PARAMETER ( MAXNLA = 999  )

      LOGICAL     OK 
      INTEGER     ILG1,ILAT,NLG1,NLAT,LONBAD,LATBAD,ICHOICE
      REAL        GLR(ILG1,ILAT),DLON(ILG1)
      REAL        GHR(NLG1,NLAT),DLAT(ILAT)

      REAL        LATM,LATP,LONM,LONP,DEGX,DEGY,TEMP,PI
      REAL        FACTLAT,FACTL,SUM,SLATP,SLATM,SLONM,SLONP
      REAL        AREA,BIGAREA,DISTRY,DISTLON,DISTMIN,DISTLAT
      REAL        RLAT(MAXNLA),RLON(MAXNLG)
      REAL        VAL(MAXNLG), DIST(MAXNLG)

      INTEGER     I,J,II,JJ,LL,IJK,N,NLG,ILG,NLPTS
      INTEGER     NEWMAX,NT,NUM,NVAL,NPP,NCHOICE
      INTEGER     LOCAT(MAXNLG), LATL(MAXNLA),LATH(MAXNLA)
      INTEGER     LOCFST(MAXNLG),IVAL(MAXNLG),NMAX(MAXNLG)
      INTEGER     LONL(MAXNLG),  LONH(MAXNLG)

*------------------------------------------------------------------------------ 
      PI = 4.0*ATAN( 1.D0 )

      LONBAD=0 
      LATBAD=0 
      OK=.TRUE. 

                              ILG = ILG1-1
      IF (MOD( ILG1,2 ).EQ.0) ILG = ILG1 
                              NLG = NLG1-1
      IF (MOD( NLG1,2 ).EQ.0) NLG = NLG1 

***    DEFINE LONGITUDE AND LATITUDE VECTORS FOR HIGH-RESOLUTION GRID. 

      DEGX=360./(FLOAT(NLG1-1)) 
      DEGY=180./(FLOAT(NLAT-1)) 

      DO  10 I=1,NLG1 
   10 RLON(I)=DEGX*FLOAT(I-1) 
      DO  20 J=1,NLAT 
   20 RLAT(J)=-90.+DEGY*FLOAT(J-1) 

***    FIND HIGH-RESOLUTION GRID POINTS WITHIN EACH LOW-RESOLUTION GRID SQUARE. 
***    THIS REQUIRES SPECIAL HANDLING OF POINTS NEAR THE POLES OR AT THE 
***    GREENWICH MERIDIAN. 

      DO 40 J=1,ILAT 
        IF(J.EQ.ILAT) THEN 
          LATP=90. 
        ELSE 
          LATP=0.5*(DLAT(J+1)+DLAT(J)) 
        ENDIF 
        IF(J.EQ.1) THEN 
          LATM=-90. 
        ELSE 
          LATM=0.5*(DLAT(J-1)+DLAT(J)) 
        ENDIF 
        DO 30 JJ=1,NLAT 
          IF(JJ.NE.1) THEN 
            IF(RLAT(JJ-1).LT.LATM.AND.RLAT(JJ).GE.LATM) LATL(J)=JJ 
          ELSE 
            IF(RLAT(JJ).GE.LATM) LATL(J)=JJ 
          ENDIF 
          IF(JJ.NE.NLAT) THEN 
            IF(RLAT(JJ).LE.LATP.AND.RLAT(JJ+1).GT.LATP) LATH(J)=JJ 
          ELSE 
            IF(RLAT(JJ).LE.LATP) LATH(J)=JJ 
          ENDIF 
   30   CONTINUE 
   40 CONTINUE 

      DO 75 I=1,ILG 
        IF(I.EQ.1) THEN 
          LONM=0.5*(DLON(ILG)+360.) 
          LONP=0.5*(DLON(I+1)+DLON(I)) 
          IF(LONM.GT.RLON(NLG)) LONM=0. 
        ELSE IF (I.EQ.ILG) THEN
          LONM=0.5*(DLON(I-1)+DLON(I)) 
          LONP=0.5*(DLON(ILG)+360.)      
        ELSE
          LONM=0.5*(DLON(I-1)+DLON(I)) 
          LONP=0.5*(DLON(I+1)+DLON(I)) 
        ENDIF 
        IF(RLON(1).GE.LONM)                           LONL(I)=1
        IF(RLON(1).LE.LONP.AND.RLON(2).GT.LONP)       LONH(I)=1
        DO 50 II=2,NLG-1
          IF(RLON(II).GE.LONM.AND.RLON(II-1).LT.LONM) LONL(I)=II 
          IF(RLON(II).LE.LONP.AND.RLON(II+1).GT.LONP) LONH(I)=II 
   50   CONTINUE 
        IF(RLON(NLG).GE.LONM.AND.RLON(NLG-1).LT.LONM) LONL(I)=NLG
        IF(RLON(NLG).LE.LONP)                         LONH(I)=NLG
        IF(LONL(I)  .GT.LONH(I)) LONH(I)=LONH(I)+NLG 
   75 CONTINUE 

***    LOOP OVER POINTS IN THE LOW-RESOLUTION FIELD, KEEPING TRACK OF ALL 
***    HIGH-RESOLUTION POINTS WITHIN EACH LOW-RESOLUTION GRID SQUARE 
***    HAVING THE SAME GROUND COVER AS THE LOW-RESOLUTION GRID SQUARE. 

      DO 900 J=1,ILAT 
        FACTL=(COS(DLAT(J)*PI/180.))**2 
        DO 800 I=1,ILG 
          NLPTS=0 
          IF(ICHOICE.EQ.1) THEN 

***          FOR CONTINUOUS FIELDS, AREA-AVERAGE THESE POINTS. 

            BIGAREA=0. 
            SUM=0. 
            DO 200 JJ=LATL(J),LATH(J) 
              IF(JJ.EQ.NLAT) THEN 
                SLATP=90. 
              ELSE 
                SLATP=0.5*(RLAT(JJ+1)+RLAT(JJ)) 
              ENDIF 
              IF(JJ.EQ.1) THEN 
                SLATM=-90. 
              ELSE 
                SLATM=0.5*(RLAT(JJ-1)+RLAT(JJ)) 
              ENDIF 
              FACTLAT=(COS(RLAT(JJ)*PI/180.))*(SLATP-SLATM)*PI/180. 
              DO 100 LL=LONL(I),LONH(I) 
                IF(LL.GT.NLG) THEN 
                  II=LL-NLG 
                ELSE 
                  II=LL 
                ENDIF 
                NLPTS=NLPTS+1 
                IF(II.EQ.1) THEN 
                  SLONM=0.5*(RLON(NLG)+360.)-360. 
                  SLONP=0.5*(RLON(II+1)+RLON(II)) 
                ELSE IF(II.EQ.NLG) THEN
                  SLONM=0.5*(RLON(II-1)+RLON(II)) 
                  SLONP=0.5*(RLON(NLG)+360.)
                ELSE 
                  SLONM=0.5*(RLON(II-1)+RLON(II)) 
                  SLONP=0.5*(RLON(II+1)+RLON(II)) 
                ENDIF 
                AREA=FACTLAT*((SLONP-SLONM)*PI/180.) 
                SUM=SUM+AREA*GHR(II,JJ) 
                BIGAREA=BIGAREA+AREA 
  100         CONTINUE 
  200       CONTINUE 
            IF(NLPTS.NE.0) THEN 
              GLR(I,J)=SUM/BIGAREA 
            ENDIF 
          ELSE 

***          FOR DISCRETE-VALUED FIELDS, KEEP TRACK OF EACH QUALIFYING HIGH- 
***          RESOLUTION POINT VALUE IN ARRAY VAL, AND THEIR DISTANCE FROM THE 
***          LOW-RESOLUTION GRID POINT IN ARRAY DIST. 

            DO 400 JJ=LATL(J),LATH(J) 
              DISTLAT=DLAT(J)-RLAT(JJ) 
ccc           print*,'jj,dlat,rlat,distlat= ',jj,dlat(j),rlat(jj),distlat
              DO 300 LL=LONL(I),LONH(I),1 
                IF(LL.GT.NLG) THEN 
                  II=LL-NLG 
                ELSE 
                  II=LL 
                ENDIF 
                NLPTS=NLPTS+1 
                VAL(NLPTS)=GHR(II,JJ) 
                DISTRY=ABS(DLON(I)-RLON(II)) 
                DISTLON=AMIN1(DISTRY,(360.-DISTRY)) 
                DIST(NLPTS)=(FACTL*(DISTLON**2)+DISTLAT**2)* 
     +                      ((PI/180.)**2) 
  300         CONTINUE 
  400       CONTINUE 
          ENDIF 

***        RETURN WITH OK=.FALSE. IF NO HIGH-RESOLUTION POINTS IN 
***        LOW-RESOLUTION GRID SQUARE DEFINED BY (LONBAD,LATBAD)=(I,J). 

          IF(NLPTS.EQ.0) THEN 
            OK=.FALSE. 
            LONBAD=I 
            LATBAD=J 
            RETURN 
          ENDIF 

***        IF ONLY ONE HIGH-RESOLUTION GRID POINT EXISTS, RETURN ITS VALUE. 

          IF(ICHOICE.NE.1.AND.NLPTS.EQ.1) THEN 
            GLR(I,J)=VAL(1) 
          ELSE IF(ICHOICE.NE.1.AND.NLPTS.GT.1) THEN 

***          CALCULATE NUMBER OF DISTINCT-VALUED POINTS (DEFINED BY NT), THEIR 
***          FIRST-OCCURRING LOCATION (DEFINED BY ARRAY LOCFST) AND THE 
***          ASSOCIATED LOCATIONS WHERE THEY OCCUR IN ARRAY VAL (DEFINED BY 
***          ARRAY NMAX). 

            NEWMAX=0 
            IF(ICHOICE.EQ.2) THEN 

***            THE FIELD BEING SMOOTHED IS GROUND COVER AND REQUIRES SPECIAL 
***            TREATMENT FOR LOW-RESOLUTION GRID SQUARES HAVING ALL THREE TYPES 
***            OF GROUND COVER. IF THE MOST FREQUENTLY-OCCURRING VALUE IS LAND 
***            BUT THERE ARE MORE (WATER+ICE) POINTS THAN LAND, THE RESULTING 
***            LOW-RESOLUTION POINT SHOULD BE RETURNED AS THE MORE FREQUENT 
***            VALUE OF WATER VERSUS ICE. 

              NT=0 
              CALL WHENEQ(NLPTS,VAL,1, 0.,IVAL,NVALW) 
              IF(NVALW.NE.0) THEN 
                NT=NT+1 
                NMAX(NT)=NVALW 
                LOCFST(NT)=IVAL(1) 
              ENDIF 
              CALL WHENEQ(NLPTS,VAL,1,+1.,IVAL,NVALI) 
              IF(NVALI.NE.0) THEN 
                NT=NT+1 
                NMAX(NT)=NVALI 
                LOCFST(NT)=IVAL(1) 
              ENDIF 
              CALL WHENEQ(NLPTS,VAL,1,-1.,IVAL,NVALL) 
              IF(NVALL.NE.0) THEN 
                NT=NT+1 
                NMAX(NT)=NVALL 
                LOCFST(NT)=IVAL(1) 
              ENDIF 
              MAXINT=MAX0(NVALL,NVALW) 
              NEWMAX=MAX0(MAXINT,NVALI) 
              IF(NT   .EQ.3             .AND.
     +          (NVALL.LT.(NVALW+NVALI)).AND.
     +           NVALL.EQ. NEWMAX) THEN 
                NT=2 
                NEWMAX=MAX0(NVALW,NVALI) 
              ELSE IF(NT   .EQ. 3            .AND.
     +               (NVALL.EQ.(NVALW+NVALI)).AND.
     +                NVALL.EQ. NEWMAX) THEN 
                NT=3 
                NMAX(1)=NVALL 
                NMAX(2)=NVALL 
                NMAX(3)=NVALL 
              ENDIF 
            ELSE 

***            DO THE GENERAL CALCULATION. 

              DO 500 N=1,NLPTS 
                IF(N.EQ.1) THEN 
                  NT=0 
                  NUM=0 
                ELSE 
                  CALL WHENEQ(N-1,VAL,1,VAL(N),IVAL,NUM) 
                ENDIF 
                IF(NUM.EQ.0) THEN 
                  NT=NT+1 
                  CALL WHENEQ(NLPTS,VAL,1,VAL(N),IVAL,NVAL) 
                  LOCFST(NT)=N 
                  NMAX(NT)=NVAL 
                  IF(NVAL.GT.NEWMAX) NEWMAX=NVAL 
                ENDIF 
  500         CONTINUE 
            ENDIF 

***          IF ONLY ONE VALUE EXISTS IN THE LOW-RESOLUTION GRID SQUARE, 
***          RETURN ITS VALUE. 

            IF(NT.EQ.1) THEN 
              GLR(I,J)=VAL(LOCFST(1)) 
            ELSE 

***            IF ONLY ONE UNIQUE-VALUED POINT HAS THE MOST FREQUENT OCCURRENCE 
***            NEWMAX, RETURN ITS VALUE. 

              CALL WHENEQI(NT,NMAX,1,NEWMAX,IVAL,NVAL) 
              IF(NVAL.EQ.1) THEN 
                GLR(I,J)=VAL(LOCFST(IVAL(1))) 
              ELSE 

***              OTHERWISE, SEARCH THROUGH THE UNIQUE-VALUED POINTS TO 
***              DETERMINE THE VALUES HAVING THE MOST FREQUENT OCCURRENCE AS 
***              NEWMAX AND THEIR LOCATION IN ARRAY VAL (DEFINED BY ARRAY 
***              LOCAT). 

                NPP=0 
                DO 600 N=1,NT 
                  IF(NMAX(N).EQ.NEWMAX) THEN 
                    TEMP=VAL(LOCFST(N)) 
                    CALL WHENEQ(NLPTS,VAL,1,TEMP,IVAL,NVAL) 
                    DO 550 IJK=1,NVAL 
                      NPP=NPP+1 
                      LOCAT(NPP)=IVAL(IJK) 
  550               CONTINUE 
                  ENDIF 
  600           CONTINUE 

***              CHOOSE THE CLOSEST POINT TO THE TARGET HAVING THAT MAXIMUM- 
***              OCCURRING VALUE. 

                DISTMIN=1.E20 
                DO 700 N=1,NPP 
                  IF(DIST(LOCAT(N)).LT.DISTMIN) THEN 
                    DISTMIN=DIST(LOCAT(N)) 
                    NCHOICE=LOCAT(N) 
                  ENDIF 
  700           CONTINUE 
                GLR(I,J)=VAL(NCHOICE) 
              ENDIF 
            ENDIF 
          ENDIF 
  800   CONTINUE 

***      REPEAT GREENWHICH MERIDIAN AS EXTRA CYCLIC LONGITUDE. 

        IF (MOD(ILG1,2).NE.0) GLR(ILG1,J) = GLR(1,J)

  900 CONTINUE 

      RETURN 
*-------------------------------------------------------------------- 

      END 
