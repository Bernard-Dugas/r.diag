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
C     $Log: pael.ftn,v $
C     Revision 3.10  2021/01/22 16:16  dugas
C     Argument PS de NIVCAL dans PAEL en 64 bits.
C
C     Revision 3.9  2016/10/26 15:25  dugas
C     Argument PSRF --> LNSP dans ELAEL (optimisation)
C
C     Revision 3.8  2014/09/25 18:42:03  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.7  2013/03/21 20:46:11  bernard
C     Corriger la documention de la routine PAEL.
C
C     Revision 3.6  2010/07/21 16:35:39  dugas
C     S'assurer que les calculs utilisent des coordonnees verticales
C      ayant la meme orientation (croissante/decroissante) quitte a
C      inverser certaines donnees a l'entree et/ou les resultats
C      a la fin (dans les routines EAPL, ELAEL, GEMAPL et PAEL).
C
C     Revision 3.5  2010/02/12 22:51:46  dugas
C     - Le code de ELAEL reflete maintenant beaucoup plus celui de PAEL.
C     - Les champs de travail sont alloues automatiquement dans les routines
C       PAEL et ELAEL, plutot que passes en argument.
C
C     Revision 3.4  2002/08/20 18:59:53  dugas
C     Modifier ELAEL pour que les coordonnees GEM d'entree et de
C        sortie puissent etre differentes.
C
C     Revision 3.3  2000/07/21 16:37:48  armnrbd
C     Ajouter la routine ELAEL (SIGMA a SIGMA).
C
C     Revision 3.2  1998/09/23 16:13:52  armnrbd
C     Ajouter les parametres TOP,BOT et INC.
C
C     Revision 3.1  1994/11/17  14:13:56  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:56:00  13:56:00  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:32:05  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.0  92/02/21  11:34:10  armnrbd
C     Initial revision
C     

      SUBROUTINE pael (FS, LA,SIG,NSL,
     +                 FP, PRLOG,NPL,
     +                 PSLOG,RLUP,RLDN, A,B )
    
      IMPLICIT none

***    FEB 02/88 - R.LAPRISE.

***    INTERPOLATES MULTI-LEVEL SET OF GRIDS FROM PRESSURE LEVELS TO
***    ETA LEVELS FOR HYBRID MODEL. INTERPOLATION IS LINEAR IN LN(PRES).

***    FP         = INPUT  GRIDS ON PRESSURE LEVELS.
***    FS         = OUTPUT GRIDS ON ETA      LEVELS.
***    PSLOG      = INPUT  GRID OF LN(SURFACE PRESSURE IN MB).
***    PRLOG(NPL) = VALUES OF INPUT PRESSURE LEVEL (Pa).
***    SIG(NSL)   = VALUES OF ETA LEVELS OF OUTPUT FIELD.
***    RLUP       = LAPSE RATE USED TO EXTRAPOLATE ABOVE PRLOG(1).
***    RLDN       = LAPSE RATE USED TO EXTRAPOLATE BELOW PRLOG(NSL).
***                 UNITS OF RLUP,RLDN ARE  DF/D(LN PRES).

***    Notes:   1) GRIDS ALL HAVE THE SAME HORIZONTAL SIZE (LA POINTS).
***             2) FP AND FS MAY BE EQUIVALENCED IN CALLING PROGRAM
***             3) PRLOG AND SIG MUST BE MONOTONIC

      INTEGER  LA,NSL,NPL,NPL1
      REAL     FP(LA,NPL),FS(LA,NSL),PSLOG (LA)
      REAL     PRLOG(NPL),RLUP,RLDN, SIG(NSL)
      REAL     A(NSL),B(NSL)
  
***    Local variables and work space.

      LOGICAL  MONOTON, REVERSE
      INTEGER  I,K,L,N, INTVL, TOPP,BOTP,INCP, TOPS,BOTS,INCS
      REAL     FPC(NPL),DFDLNP(0:NPL+1),DLNP(NPL),X, LOCSIG(NSL)
      REAL     FSC(NSL),HOLD

      LOGICAL           INFO
      COMMON  /ZZVERBO/ INFO

      EXTERNAL MONVERT,NIVCAL,XIT

*---------------------------------------------------------------------

      LOCSIG(:) = SIG(:)

***    CHECK THAT PRLOG IS MONOTONIC.

      CALL MONVERT( PRLOG,NPL, TOPP,BOTP,INCP, MONOTON )

      IF (.NOT.MONOTON)                                        THEN
          IF (INFO) WRITE(6,6000) 'PRESSURE'
          CALL                                     XIT('  Pael  ',-1 )
      END IF

***    PRECOMPUTE INVERSE OF DELTA LN PRES.

      DO  L=TOPP,BOTP-INCP,INCP
          DLNP(L) = 1.0 / (PRLOG(L+INCP)-PRLOG(L))
      END DO
  
***    LOOP OVER ALL HORIZONTAL POINTS.
  
      DO  300 I=1,LA 
  
***        GET A COLUMN OF FP (ON PRESSURE LEVELS).
  
          FPC(:) = FP(I,:)
  
***        COMPUTE VERTICAL DERIVATIVE OVER ALL PRESSURE INTERVALS.
  
          DFDLNP(TOPP     ) = RLUP
          DO  L=TOPP,BOTP-INCP,INCP
          DFDLNP(L   +INCP) = (FPC(L+INCP)-FPC(L)) * DLNP(L)
          END DO
          DFDLNP(BOTP+INCP) = RLDN
  
***        COMPUTE THE LOCAL SIGMA VALUES. 
  
          CALL NIVCAL( LOCSIG, A,B, 100.D0 * EXP( PSLOG(I) ), NSL,1,1 )

          IF (I == 1)                                          THEN

***            CHECK THAT THE LOCAL SIGMA VECTOR IS MONOTONIC
***            (AND IN WHAT DIRECTION). DO THIS ONLY ONCE.

              CALL MONVERT( LOCSIG,NSL, TOPS,BOTS,INCS, MONOTON )

              IF (.NOT.MONOTON)                                THEN
                  IF (INFO) WRITE(6,6000) 'LOCAL SIGMA'
                  CALL                             XIT('  Pael  ',-2 )
              END IF

***            ARE THE TWO COORDINATES IN THE SAME ORDER ?

              REVERSE = .NOT. ( (TOPS < BOTS) .EQV. (TOPP < BOTP) )

          END IF

          IF (REVERSE)                                         THEN
              DO  L=1,NSL/2
                  HOLD            = LOCSIG(    L  )
                  LOCSIG(    L  ) = LOCSIG(NSL-L+1)
                  LOCSIG(NSL-L+1) = HOLD
              END DO
              IF (I == 1)                                      THEN
                  INCS =-INCS
                  L    = BOTS
                  BOTS = TOPS
                  TOPS = L
              END IF
          END IF
  
          LOCSIG(:) = ALOG( LOCSIG(:) )

***        LOOP OVER SIGMA LEVELS TO BE INTERPOLATED.
***        X IS THE LN(PRES)=LN(PS)*LN(SIG) VALUE OF REQUIRED SIGMA LEVEL.
  
          K = TOPP
          DO  200 N=TOPS,BOTS,INCS

              X = LOCSIG(N)+PSLOG(I)
  
***            FIND WHICH PRESSURE INTERVAL WE ARE IN.
  
              DO  L=K,BOTP,INCP
                  INTVL = L
                  IF (X.LT.PRLOG(L)) GOTO 100
              END DO

              INTVL = BOTP +INCP
  100         K     = INTVL-INCP

              IF (K.EQ.TOPP-INCP) K = TOPP
  
***            NOW INTERPOLATE AT THIS POINT.
  
              FSC(N)  = FPC(K)+DFDLNP(INTVL)*(X-PRLOG(K))

  200     CONTINUE

          IF (REVERSE)                                         THEN
              DO  L=1,NSL/2
                  HOLD         = FSC(    L  )
                  FSC(    L  ) = FSC(NSL-L+1)
                  FSC(NSL-L+1) = HOLD
              END DO
          END IF

          FS(I,:) = FSC(:)
  
  300 CONTINUE
  
      RETURN
*---------------------------------------------------------------------

 6000 FORMAT(' Pael: ',A,' levels are not monotonic.'/)

      END 
      SUBROUTINE elael( FP, PR ,NPL, AI,BI, COORDI,PTOITI,
     +                  FS, SIG,NSL, AO,BO, COORDO,PTOITO, 
     +                  LNSP,PTOP,PTOPI, RLUP,RLDN, LA )
    
      IMPLICIT    none

***    JUNE 2000 - NILS EK
***    FEB 02/88 - R.LAPRISE.

***    INTERPOLATES MULTI-LEVEL SET OF GRIDS BETWEEN DIFFERENT SETS OF
***    ETA LEVELS. INTERPOLATION IS LINEAR IN LN(PRES).

***    FP         = INPUT  GRIDS ON ETA LEVELS.
***    FS         = OUTPUT GRIDS ON ETA LEVELS.
***    LNSP       = INPUT  GRID OF SURFACE LN PRESSURE IN Pa
***    PR(NPL)    = VALUES OF INPUT ETA LEVEL
***    SIG(NSL)   = VALUES OF ETA LEVELS OF OUTPUT FIELD.
***    PTOP       = PRESSURE GRID AT LID IN PA (USED ONLY FOR GEM OUTPUT)
***    PTOPI      = PRESSURE GRID AT LID IN PA (USED ONLY FOR GEM INPUT)
***    RLUP       = LAPSE RATE USED TO EXTRAPOLATE ABOVE PRLOG(1).
***    RLDN       = LAPSE RATE USED TO EXTRAPOLATE BELOW PRLOG(NSL).
***                 UNITS OF RLUP,RLDN ARE  DF/D(LN PRES).

***    Notes:   1) GRIDS ALL HAVE THE SAME HORIZONTAL SIZE (LA POINTS).
***             2) PR AND SIG MUST BE MONOTONIC

      CHARACTER*4 COORDI,COORDO
      INTEGER     LA,NSL,NPL,NPL1
      REAL        PR (NPL),AI(NPL),BI(NPL),PTOITI
      REAL        SIG(NSL),AO(NSL),BO(NSL),PTOITO 
      REAL        FS(LA,NSL),PTOP(LA),PTOPI(LA)
      REAL        FP(LA,NPL),LNSP(LA)
      REAL        RLUP,RLDN
  
***    Local variables

      REAL(8)     CORR2
      LOGICAL     MONOTON,REVERSE
      REAL(8)     ZWB,ZWT,DFDLNP(0:NPL+1),FSC(NSL)
      REAL(8)     FPC(NPL),PI(NPL),PO(NSL),HOLD
      INTEGER     I,KI,KO,L, TOPP,BOTP,INCP
      INTEGER     INTVL, TOPS,BOTS,INCS

      REAL(8),    DIMENSION(:), POINTER, SAVE :: PSRF

      LOGICAL              INFO
      COMMON     /ZZVERBO/ INFO

      EXTERNAL    MONVERT,XIT

*---------------------------------------------------------------------
      CORR2 = LOG( 100000._8 )

      IF ((COORDI /= 'GEM4') .OR. (COORDO /= 'GEM4'))          THEN
          IF (.NOT.ASSOCIATED( PSRF )) ALLOCATE( PSRF(LA) )
          PSRF = EXP( 1.0_8 * LNSP )
      END IF

***    CHECK THAT PR IS MONOTONIC.

      CALL MONVERT( PR,NPL, TOPP,BOTP,INCP, MONOTON )

      IF (.NOT.MONOTON)                                        THEN
          IF (INFO) WRITE(6,6000) COORDI
          CALL                                     XIT(' Elael  ',-1 )
      END IF
 
***    CHECK THAT SIG IS MONOTONIC.

      CALL MONVERT( SIG,NSL, TOPS,BOTS,INCS, MONOTON )

      IF (.NOT.MONOTON)                                        THEN
          IF (INFO) WRITE(6,6000) COORDO
          CALL                                     XIT(' Elael  ',-2 )
      END IF
 
***    ARE THE TWO COORDINATES IN THE SAME ORDER ?

      REVERSE = .NOT. ( (TOPS < BOTS) .EQV. (TOPP < BOTP) )

      IF (REVERSE)                                             THEN
          INCS =-INCS
          L    = BOTS
          BOTS = TOPS
          TOPS = L
      END IF

***    LOOP OVER ALL HORIZONTAL POINTS.
  
      DO  300 I=1,LA 
  
***        GET A COLUMN OF FP (ON ORIGINAL ETA LEVELS).
  
          FPC(:) = FP(I,:)
  
***        COMPUTE LOCAL PRESSURE AT EACH POINT FOR BOTH SETS OF ETA
***        FIRST, THE INPUT LEVELS

          IF (COORDI.EQ.'GEM')                                 THEN

              DO  L=1,NPL
                  PI(L) = PR(L) * (PSRF(I) - PTOPI(I) ) + PTOPI(I)
              END DO

          ELSE IF (COORDI /= 'GEM4')                           THEN

              DO  L=1,NPL
                  PI(L) = AI(L) + BI(L) * PSRF(I)
              END DO

          ELSE

              DO  L=1,NPL ! Calculate LN( local pressure )
                  PI(L) = AI(L) + BI(L) * (LNSP(I) - CORR2)
              END DO

          END IF

***        COMPUTE INPUT VERTICAL DERIVATIVE OVER ALL PRESSURE INTERVALS
***        AND CONVERT PI FROM P TO LOG(P).
  
          DFDLNP(TOPP)      = RLUP
          DFDLNP(BOTP+INCP) = RLDN

          IF (COORDI /= 'GEM4') PI(:) = LOG( PI(:) )

          DO  L=TOPP,BOTP-INCP,INCP
              DFDLNP(L+INCP) = (FPC(L+INCP)-FPC(L)) / (PI(L+INCP)-PI(L))
          END DO
  
***        NEXT, THE OUTPUT LEVELS

          IF (COORDO.EQ.'GEM')                                 THEN

              DO  L=1,NSL
                  PO(L) = SIG(L) * ( PSRF(I) - PTOP(I) ) + PTOP(I)
              END DO

          ELSE IF (COORDO /= 'GEM4')                           THEN

              DO  L=1,NSL
                  PO(L) = AO(L) + BO(L) * PSRF(I)
              END DO

          ELSE

              DO  L=1,NSL ! Calculate LN( local pressure )
                  PO(L) = AO(L) + BO(L) * (LNSP(I) - CORR2)
              END DO

          END IF

          IF (REVERSE)                                         THEN
              DO  L=1,NSL/2
                  HOLD        = PO(    L  )
                  PO(    L  ) = PO(NSL-L+1)
                  PO(NSL-L+1) = HOLD
              END DO
          END IF
  
          IF (COORDO /= 'GEM4') PO(:) = LOG( PO(:) )

***        PERFORM VERTICAL INTERPOLATION OF FP -> FS BETWEEN THE TWO 
***        SETS OF PRESSURE LEVELS.

          KI = TOPP 

          DO  KO = TOPS,BOTS,INCS

***            FIND WHICH PRESSURE INTERVAL WE ARE IN.
  
              DO  L=KI,BOTP,INCP
                  INTVL = L
                  IF (PO(KO) < PI(L)) GOTO 100
              END DO

              INTVL = BOTP +INCP
  100         KI    = INTVL-INCP

              IF (KI == TOPP-INCP) KI = TOPP
  
***            NOW INTERPOLATE AT THIS POINT.

              FSC(KO) = FPC(KI)+DFDLNP(INTVL)*(PO(KO)-PI(KI))

          END DO

          IF (REVERSE)                                         THEN
              DO  L=1,NSL/2
                  HOLD         = FSC(    L  )
                  FSC(    L  ) = FSC(NSL-L+1)
                  FSC(NSL-L+1) = HOLD
              END DO
          END IF

          FS(I,:) = FSC(:)
  
  300 CONTINUE
  
      RETURN
*---------------------------------------------------------------------

 6000 FORMAT(' Elael: ',A,' levels are not monotonic.'/)

      END 

