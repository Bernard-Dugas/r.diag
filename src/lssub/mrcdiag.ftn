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
*     $Log: mrcdiag.ftn,v $
*     Revision 3.5  2014/09/25 18:42:03  dugas
*     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
*
*     Revision 3.4  1998/09/18 01:41:12  armnrbd
*     Modifier GAPL pour qu'elle fonctionne avec des
*     coordonnees verticales croissantes ou decroissantes.
*
*     Revision 3.3  1995/12/06  20:40:25  armnrbd
*     Ajouter la routine SDET.
*
*     Revision 3.2  1995/11/22  02:27:59  armnrbd
*     Ajouter la routine CLATLON.
*
*     Revision 3.1  1995/11/08  10:19:12  armnrbd
*     Ajouter la routine DEFCPHY.
*
*     Revision 3.0  95/11/07  01:22:37  01:22:37  armnrbd (Bernard Dugas)
*     Version initiale.
*     

      SUBROUTINE BX0 (F,G,FNI,FNJ,GNI,GNJ,INTERP)

      IMPLICIT none

***    NOV 04/93 - M. GIGUERE  (INTERPOLATION CUBIQUE.)
***    JUL 08/93 - B.DENIS, S.TURNER, S.BINER

***    CE SOUS-PROGRAMME FAIT LA MOYENNE ENTRE DEUX POINTS DU TABLEAU F 
***    CONSECUTIFS SELON L'AXE DES X, ET MET CETTE VALEUR DANS LE TABLEAU G.
***    SI INTERP=4 ALORS INTERPOLATION CUBIQUE EN X.

      INTEGER  FNI,FNJ,GNI,GNJ,INTERP
      REAL     F(FNI,FNJ), G(GNI,GNJ)

      REAL     CF
      INTEGER  I,J,ISB

*---------------------------------------------------
      DO  J=1,GNJ
          DO  I=1,GNI
              G(I,J) = 0.5 * ( F(I,J)+F(I+1,J) )
          END DO
      END DO

      IF (INTERP.NE.4) RETURN
      ISB = 1
      CF  = -.125

***    APPEL A UN OPERATEUR A TROIS POINTS
***    TYPE: FX= FX + CF * [ F(X+DX)+F(X-DX)-2*FX ]

      CALL XFILT( G,G,GNI,GNJ,1,ISB,0,0,CF )

      RETURN
*---------------------------------------------------

      END
      SUBROUTINE BX1 (F,G,FNI,FNJ,GNI,GNJ,INTERP)

      IMPLICIT none

***    NOV 04/93 - M. GIGUERE (INTERPOLATION CUBIQUE)
***    JUL 08/93 - S.TURNER, S.BINER, B.DENIS

***    CE SOUS-PROGRAMME FAIT LA MOYENNE ENTRE DEUX POINTS DU TABLEAU F 
***    CONSECUTIFS SELON L'AXE DES X, ET MET CETTE VALEUR DANS LE TABLEAU G.
***    DE PLUS, LES EXTREMITES SONT COMBLEES EN Y METTANT LA VALEUR LA PLUS 
***    PROCHE LE LONG DE L'AXE DES X.
***    SI INTERP=4 ALORS INTERPOLATION CUBIQUE EN X.

      INTEGER  FNI,FNJ,GNI,GNJ,INTERP
      REAL     F(FNI,FNJ), G(GNI,GNJ)

      REAL     CF
      INTEGER  I,J,ISB,IS,JS,KS

*---------------------------------------------------
      DO  J=1,GNJ
          DO  I=1,GNI-2
              G(I+1,J) = 0.5 * ( F(I,J)+F(I+1,J) )
          END DO
      END DO

      IF (INTERP.NE.4) GOTO 99

      IS  = 2
      JS  = 0
      KS  = 0
      ISB = IS+1
      CF  = -.125

***    APPEL A UN OPERATEUR A TROIS POINTS
***    TYPE: FX= FX + CF * [ F(X+DX)+F(X-DX)-2*FX ]

      CALL XFILT( G(2,1),G(2,1),GNI,GNJ,1,ISB,JS,KS,CF )

 99   CONTINUE

      DO  J=1,GNJ

***        REMPLISSAGE DES FRONTIERES

          G(1,  J) = G(2,    J)
          G(GNI,J) = G(GNI-1,J)

      END DO

      RETURN
*---------------------------------------------------

      END
      SUBROUTINE BY0 (F,G,FNI,FNJ,GNI,GNJ,INTERP)

      IMPLICIT none

***    NOV 04/93 - M. GIGUERE  (INTERPOLATION CUBIQUE.)
***    JUL 08/93 - S.TURNER, B.DENIS, S.BINER

***    CE SOUS-PROGRAMME FAIT LA MOYENNE ENTRE DEUX POINTS DU TABLEAU F 
***    CONSECUTIFS SELON L'AXE DES Y, ET MET CETTE VALEUR DANS LE TABLEAU G.
***    SI INTERP=4 ALORS INTERPOLATION CUBIQUE EN Y.

      INTEGER  FNI,FNJ,GNI,GNJ,INTERP
      REAL     F(FNI,FNJ), G(GNI,GNJ)

      REAL     CF
      INTEGER  I,J,JSB

*---------------------------------------------------
      DO  I=1,GNI
          DO  J=1,GNJ
              G(I,J) = 0.5 * ( F(I,J)+F(I,J+1) ) 
          END DO
      END DO

      IF (INTERP.NE.4) RETURN
      JSB = 1
      CF  = -.125

***    APPEL A UN OPERATEUR A TROIS POINTS
***    TYPE: FY= FY + CF * [ F(Y+DY)+F(Y-DY)-2*FY ]

      CALL YFILT( G,G,GNI,GNJ,1,0,JSB,0,CF )

      RETURN
*---------------------------------------------------

      END
      SUBROUTINE BY1 (F,G,FNI,FNJ,GNI,GNJ,INTERP)

      IMPLICIT none

***    NOV 04/93 - M. GIGUERE (INTERPOLATION CUBIQUE)
***    JUL 08/93 - B.DENIS, S.TURNER, S.BINER

***    CE SOUS-PROGRAMME FAIT LA MOYENNE ENTRE DEUX POINTS DU TABLEAU F 
***    CONSECUTIFS SELON L'AXE DES Y, ET MET CETTE VALEUR DANS LE TABLEAU G.
***    DE PLUS, LES EXTREMITES SONT COMBLEES EN Y METTANT LA VALEUR LA PLUS
***    PROCHE LE LONG DE L'AXE DES Y.
***    SI INTERP=4 ALORS INTERPOLATION CUBIQUE EN Y.

      INTEGER  FNI,FNJ,GNI,GNJ,INTERP
      REAL     F(FNI,FNJ), G(GNI,GNJ)

      REAL     CF
      INTEGER  I,J,JSB,IS,JS,KS

*---------------------------------------------------
      DO  I=1,GNI
          DO  J=1,GNJ-2
              G(I,J+1) = 0.5 * ( F(I,J)+F(I,J+1) )
          END DO
      END DO

      IF (INTERP.NE.4) GOTO 99

      IS  = 0
      JS  = 2
      KS  = 0
      JSB = JS+1
      CF  = -.125

***    APPEL A UN OPERATEUR A TROIS POINTS
***    TYPE: FY= FY + CF * [ F(Y+DY)+F(Y-DY)-2*FY ]

      CALL YFILT( G(1,2),G(1,2),GNI,GNJ,1,IS,JSB,KS,CF )

 99   CONTINUE

      DO  I=1,GNI

***        REMPLISSAGE DES FRONTIERES

         G(I,  1) = G(I,    2)
         G(I,GNJ) = G(I,GNJ-1)

      END DO

      RETURN
*---------------------------------------------------

      END
      SUBROUTINE CALLNPH (LNPHT,T,DZHTM,LNPM,NIM,NJM,NK,NKP,GRAV,RGAS)

      IMPLICIT none

***    M.GIGUERE                            29 NOVEMBRE 1993

***    CE S-P CALCULE LE LOGARITHME NATUREL DE LA PRESSION
***    HYDROSTATIQUE AU TOIT
***    HYPOTHESE: TEMPERATURE ISOTHERME AU DESSUS DU DERNIER
***    NIVEAU MOMENTUM.

      INTEGER  NIM,NJM,NK,NKP
      REAL     DZHTM(NIM,NJM),LNPM(NIM,NJM,NKP)
      REAL     T(NIM,NJM,NK),LNPHT(NIM,NJM)
      REAL     GRAV,RGAS

      INTEGER  I,J

*-----------------------------------------------------------------------
      DO  J=1,NJM
          DO  I=1,NIM
              LNPHT(I,J) = LNPM(I,J,1)-(GRAV/(RGAS*T(I,J,1)))*DZHTM(I,J)
          END DO
      END DO

      RETURN
*-----------------------------------------------------------------------

      END
      SUBROUTINE CALLNPM (T,DZPM,LNSP,LNPM,NIM,NJM,NK,NKP,GRAV,RGAS)

      IMPLICIT none

***    MAR 22/94 - M.GIGUERE (COSMETIQUES)
***    JUL  1/93 - G.BERGERON et C.THURRE

***    CE S-P CALCULE LE LOGARITHME NATUREL DE LA PRESSION
***    HYDROSTATIQUE AUX POINTS DE GRILLE SITUES SUR DES
***    NIVEAUX DE TYPE MOMENTUM.

      INTEGER  NIM,NJM,NK,NKP
      REAL     DZPM(NIM,NJM,NK),LNSP(NIM,NJM),LNPM(NIM,NJM,NKP)
      REAL     T(NIM,NJM,NK),GRAV,RGAS

      INTEGER  I,J,K

*-----------------------------------------------------------------------
***    LES CONDITIONS LIMITES POUR 'LNPM' SONT DONNEES PAR 
***    'LNSP', I.E. LE LOG. NAT. DE LA PRESION DE SURFACE.

      DO  J=1,NJM
          DO  I=1,NIM
              LNPM(I,J,NKP) = LNSP(I,J)
          END DO
      END DO

***    CALCULONS 'LNPM'

      DO  K=NK,1,-1
          DO  J=1,NJM
              DO  I=1,NIM
                  LNPM(I,J,K) = LNPM(I,J,K+1)
     +                        - GRAV*DZPM(I,J,K) / (RGAS*T(I,J,K))
              END DO
          END DO
      END DO

      RETURN
*-----------------------------------------------------------------------

      END
      SUBROUTINE CALLNPT (LNPT,LNPM,T,ZPT,ZPM,NIM,NJM,NK,NKP,
     +                    GRAV,RGAS,HTOIT)

      IMPLICIT none

***    DEC  1/93 - M.GIGUERE (CALCUL EXACT DE LA PRESSION)
***    JUL  1/93 - G.BERGERON ET C.THURRE

***    CE S-P INTERPOLE LE LOGARITHME NATUREL DE LA PRESSION
***    HYDROSTATIQUE AUX POINTS DE GRILLE SITUES SUR DES
***    NIVEAUX DE TYPE THERMODYNAMIQUE.

      INTEGER  NIM,NJM,NK,NKP
      REAL     LNPM(NIM,NJM,NKP),LNPT(NIM,NJM,NK) 
      REAL     T(NIM,NJM,NK),ZPT(NIM,NJM,NK),ZPM(NIM,NJM,NKP)
      REAL     GRAV,RGAS,HTOIT

      REAL     ROG
      INTEGER  I,J,K

*-----------------------------------------------------------------------
***    CALCULONS 'LNPT'

      ROG = GRAV/RGAS

      DO  K=1,NK
          DO  J=1,NJM
              DO  I=1,NIM
                  LNPT(I,J,K) = LNPM(I,J,K)
     +                        - (ROG/T(I,J,K))*(ZPT(I,J,K)-ZPM(I,J,K))
              END DO
          END DO
      END DO

      RETURN
*-----------------------------------------------------------------------

      END
      SUBROUTINE CALZ (H0,ZGT,ZGM,ZPT,ZPM,DZPM,NIM,NJM,NK,NKM,NKP,HTOIT)

      IMPLICIT none

***    AUG 04/94 - M.GIGUERE (CORRECTION POUR TENIR COMPTE DU PREMIER
***                           NIVEAU THERMO. DECALE OU NON.)
***    DEC  1/93 - M.GIGUERE (CORRECTION DU CALCUL DES NIVEAUX.)
***    JUL  1/93 - G.BERGERON et C.THURRE

***    CE S-P FAIT DANS L'ORDRE :
***    1) LA CALCUL DE LA HAUTEUR PHYSIQUE 'ZPT' ASSOCIEE A CHACUN DES
***       POINTS DE GRILLE POUR LES NIVEAUX GAL-CHEN DE TYPE THERMO.
***    2) L'INTERPOLATION DE LA HAUTEUR PHYSIQUE 'ZPM' SUR DES NIVEAUX
***       GAL-CHEN DE TYPE MOMENTUM.
***    3) LA DIFFERENCE D'HAUTEUR PHYSIQUE 'DZPM' ENTRE DEUX NIVEAUX
***       GAL-CHEN DE TYPE MOMENTUM.

      INTEGER  NIM,NJM,NK,NKM,NKP
      REAL     H0(NIM,NJM),ZGT(NK),ZGM(NK),ZPT(NIM,NJM,NK)
      REAL     ZPM(NIM,NJM,NKP),DZPM(NIM,NJM,NK),HTOIT

      INTEGER  I,J,K

*-----------------------------------------------------------------------
***    CALCULONS 'ZGM'

      DO  K=2,NKM
          ZGM(K) = (ZGT(K)+ZGT(K-1))*0.5
      END DO

***    CONDITIONS LIMITES POUR 'ZGM'

      ZGM(1)  = (HTOIT+ZGT(1))*0.5
      ZGM(NK) = ZGT(NK-1)/2.0

***    CALCULONS 'ZPT' ET 'ZPM'

      DO  K=1,NK
          DO  J=1,NJM
              DO  I=1,NIM
                  ZPT(I,J,K) = H0(I,J)+(1.-H0(I,J)/HTOIT)*ZGT(K)
                  ZPM(I,J,K) = H0(I,J)+(1.-H0(I,J)/HTOIT)*ZGM(K)
              END DO
          END DO
      END DO

      DO  J=1,NJM
          DO  I=1,NIM
              ZPM(I,J,NK+1) = H0(I,J)
          END DO
      END DO

***    CALCULONS 'DZPM'

      DO  K=1,NK
          DO  J=1,NJM
              DO  I=1,NIM
                  DZPM(I,J,K) = ZPM(I,J,K)-ZPM(I,J,K+1)
              END DO
          END DO
      END DO

      RETURN
*-----------------------------------------------------------------------

      END
      SUBROUTINE CLATLON (DEGLAT,DEGLON,J,GRPI,GRPJ,GRD60,
     +                    GRDGRW,KHEM,IL1,IL2,ILG)

      IMPLICIT  none

***    D.CAYA 6 AVRIL 1993

***    CALCULE LES LATITUDES ET LES LONGITUDES DES
***    POINTS DE LA J-IEME TRANCHE DE LA GRILLE PS.

      INTEGER  IL1,IL2,ILG,J,KHEM
      REAL     GRPI,GRPJ,GRD60,GRDGRW
      REAL     DEGLAT(ILG),DEGLON(ILG)

      INTEGER  I
      REAL     YPS,XPS

      EXTERNAL LLFXY

C---------------------------------------------------------------------
      YPS = FLOAT( J ) - GRPJ
      DO  I = IL1,IL2
          XPS = FLOAT( I ) - GRPI
          CALL LLFXY( DEGLAT(I),DEGLON(I),XPS,YPS,GRD60,GRDGRW,KHEM )
          IF (DEGLON(I).LT.0.) DEGLON(I) = DEGLON(I) + 360.
      END DO
C
      RETURN
C---------------------------------------------------------------------

      END
      SUBROUTINE DEFCPHY (GRAV,RGAS)

***    JUL  1/93 - G.BERGERON et C.THURRE

***    CE S-P DEFINIT LES CONSTANTES SUIVANTES:
***    GRAV = ACCELERATION GRAVITATIONNELLE (M/SEC**2).
***    RGAS = CONSTANTE DES GAZ POUR L'AIR SEC (JOULE/(KG*DEG)).

      REAL GRAV,RGAS

*-----------------------------------------------------------------------
      GRAV = 9.80616
      RGAS = 287.04

      RETURN
*-----------------------------------------------------------------------

      END
      SUBROUTINE DX0 (F,D60,G,FNI,FNJ,GNI,GNJ)

      IMPLICIT none

***    JUL 08/93 - S.TURNER, B.DENIS, S.BINER

***    CE SOUS-PROGRAMME FAIT LA DERIVEE ENTRE DEUX POINT DU TABLEAU F 
***    CONSECUTIFS SELON L'AXE DES X, ET MET CETTE DERIVEE DANS LE TABLEAU
***    G. LA DERIVEE EST FAITE EN DIVISANT LA DIFFERENCE ENTRE DEUX POINTS
***    PAR LA LONGUEUR DE REFERENCE D60.

***    EX: DF(I,J)/DX = (F(I+1,J)-F(I,J))/D60

      INTEGER  FNI,FNJ,GNI,GNJ
      REAL     F(FNI,FNJ), G(GNI,GNJ), D60

      INTEGER  I,J

*-------------------------------------------------
      DO  J=1,GNJ
          DO  I=1,GNI
              G(I,J) = ( F(I+1,J)-F(I,J)) * (1.0/D60)
          END DO
      END DO

      RETURN
*-------------------------------------------------

      END
      SUBROUTINE DX1 (F,D60,G,FNI,FNJ,GNI,GNJ)

      IMPLICIT none

***    JUL 08/93 - B.DENIS, S.BINER, S.TURNER

***    CE SOUS-PROGRAMME FAIT LA DERIVEE ENTRE DEUX POINTS DU TABLEAU F 
***    CONSECUTIFS SELON L'AXE DES X, ET MET CETTE DERIVEE DANS LE TABLEAU
***    G. DE PLUS, LES EXTREMITES SONT COMBLEES EN Y METTANT LA VALEUR LA 
***    PLUS PROCHE LE LONG DE L'AXE DES X. LA DERIVEE EST FAITE EN DIVISANT 
***    LA DIFFERENCE ENTRE DEUX POINTS PAR LA LONGUEUR DE REFERENCE D60.

***    EX: DF(I,J)/DX = (F(I+1,J)-F(I,J))/D60

      INTEGER  FNI,FNJ,GNI,GNJ
      REAL     F(FNI,FNJ), G(GNI,GNJ),D60

      INTEGER  I,J

*------------------------------------------------------
      DO  J=1,GNJ

          DO  I=1,GNI-2
              G(I+1,J) = (F(I+1,J)-F(I,J)) * (1.0/D60) 
          END DO

***        REMPLISSAGE DES FRONTIERES

          G(1,  J) = G(2,    J)
          G(GNI,J) = G(GNI-1,J)

      END DO

      RETURN
*------------------------------------------------------

      END
      SUBROUTINE DY0 (F,D60,G,FNI,FNJ,GNI,GNJ)

      IMPLICIT none

***    JUL 08/93 - S.TURNER, B.DENIS, S.BINER

***    CE SOUS-PROGRAMME FAIT LA DERIVEE ENTRE DEUX POINT DU TABLEAU F 
***    CONSECUTIFS SELON L'AXE DES Y, ET MET CETTE DERIVEE DANS LE TABLEAU
***    G. LA DERIVEE EST FAITE EN DIVISANT LA DIFFERENCE ENTRE DEUX POINTS
***    PAR LA LONGUEUR DE REFERENCE D60.

***    EX: DF(I,J)/DY = (F(I,J+1)-F(I,J))/D60

      INTEGER  FNI,FNJ,GNI,GNJ
      REAL     F(FNI,FNJ), G(GNI,GNJ), D60

      INTEGER  I,J

*--------------------------------------------------
      DO  I=1,GNI
          DO  J=1,GNJ
              G(I,J) = ( F(I,J+1)-F(I,J)) * (1.0/D60)
          END DO
      END DO

      RETURN
*--------------------------------------------------

      END
      SUBROUTINE DY1 (F,D60,G,FNI,FNJ,GNI,GNJ)

      IMPLICIT none

***    JUL 08/93 - B.DENIS, S.TURNER, S.BINER

***    CE SOUS-PROGRAMME FAIT LA DERIVEE ENTRE DEUX POINTS DU TABLEAU F
***    CONSECUTIFS SELON L'AXE DES Y, ET MET CETTE DERIVEE DANS LE TABLEAU
***    G. DE PLUS, LES EXTREMITES SONT COMBLEES EN Y METTANT LA VALEUR LA 
***    PLUS PROCHE LE LONG DE L'AXE DES Y. LA DERIVEE EST FAITE EN DIVISANT 
***    LA DIFFERENCE ENTRE DEUX POINTS PAR LA LONGUEUR DE REFERENCE D60.

***    EX: DF(I,J)/DY = (F(I,J+1)-F(I,J))/D60

      INTEGER  FNI,FNJ,GNI,GNJ
      REAL     F(FNI,FNJ), G(GNI,GNJ),D60

      INTEGER  I,J

*------------------------------------------------
      DO  I=1,GNI

          DO  J=1,GNJ-2
              G(I,J+1)= (F(I,J+1)-F(I,J)) * (1.0/D60)
          END DO

***        REMPLISSAGE DES FRONTIERES

          G(I,  1) = G(I,    2)
          G(I,GNJ) = G(I,GNJ-1)

      END DO

      RETURN
*------------------------------------------------

      END
      SUBROUTINE DZHTMOM (DZHTM,ZPM,NIM,NJM,NKP,HTOIT)

      IMPLICIT none

***    MICHEL GIGUERE                      DEC 1993

***    CALCUL L'EPAISSEUR ENTRE LE TOIT ET LE NIVEAU
***    MOMENTUM LE PLUS ELEVE.

      INTEGER  NIM,NJM,NKP,I,J
      REAL     ZPM(NIM,NJM,NKP),DZHTM(NIM,NJM)
      REAL     HTOIT

*-----------------------------------------------------------------------
      DO  I=1,NIM
          DO  J=1,NJM

              DZHTM(I,J) = HTOIT-ZPM(I,J,1)

          END DO
      END DO

      RETURN
*-----------------------------------------------------------------------

      END
      SUBROUTINE GALCPHI (PHI, T,PSLN,LNPM,
     +                   RGAS, LEN, NSL,NSLP, SIG) 
  
      IMPLICIT none

***    FEB 23/93 - M.GIGUERE.

***    COMPUTE PHI FROM T ON GAL-CHEN COORDINATES.

***    PHI ANF T MAY BE EQUIVALENCED IN THE CALLING PROGRAM. 
***    PHIS, THE SURFACE GEOPOTENTIAL, MUST BE IN PHI( ,NSL+1) 
***    PSLN IS LN(PS), FOR PS IN PA. 
***    NSLP = NSL + 1. 

      INTEGER  LEN,NSL,NSLP
      REAL     PHI(LEN,NSLP),T(LEN,NSL),LNPM(LEN,NSL)
      REAL     PSLN(LEN),SIG(NSLP),RGAS

      INTEGER  L,N

*-----------------------------------------------------------------------
      SIG(NSL+1) = 0. 
  
      DO  N=1,LEN

          DO  L=1,NSL
              SIG(L)   = LNPM(N,L)-PSLN(N)
          END DO

          DO  L=1,NSL 
              SIG(L)   = SIG(L+1)-SIG(L)
          END DO
  
          DO  L=NSL,1,-1
              PHI(N,L) = PHI(N,L+1) +RGAS*T(N,L)*SIG(L) 
          END DO
  
      END DO

      RETURN
*-----------------------------------------------------------------------

      END 
      SUBROUTINE GAPL  (FP, LA,PRLOG,NPL, FS,PRES,NSL, PSLOG,RLUP,RLDN,
     +                  NSL1,FSC,DFDLNS,DLNS) 
  
      IMPLICIT none

***    AUG 18/93 - M.GIGUERE.

***    INTERPOLATES MULTI-LEVEL SET OF GRIDS FROM GAL-CHEN LEVELS 
***    TO PRESSURE LEVELS FOR REGIONAL CLIMATE MODEL.
***    THE INTERPOLATION IS LINEAR IN LN(PRES).

***    ALL GRIDS HAVE THE SAME HORIZONTAL SIZE (LA POINTS).

***    FS          =  INPUT GRIDS ON GAL-CHEN LEVELS.
***    FP          = OUTPUT GRIDS ON PRESSURE LEVELS.
***   (NOTE THAT FP AND FS MAY BE EQUIVALENCED IN CALLING PROGRAM)

***    PSLOG       = INPUT  GRID OF LN(SURFACE PRESSURE IN PA).
***    PRLOG(NPL)  = LOG PRESSURE VALUES OF OUTPUT LEVELS (PA). 
***    PRES(LA,NSL)= LOG PRESSURE VALUES OF INPUT FIELD. 
***   (BOTH MUST BE MONOTONIC).

***    NSL1        = NSL+1.
***    RLUP        = LAPSE RATE USED TO EXTRAPOLATE ABOVE TOP ETA. 
***    RLDN        = LAPSE RATE USED TO EXTRAPOLATE BELOW BOTTOM ETA.
***                  UNITS OF RLUP,RLDN ARE  DF/D(LN PRES).

      INTEGER  LA,NPL,NSL,NSL1
      REAL     FP(LA,NPL), FS(LA,NSL), PSLOG(LA) 
      REAL     PRLOG(NPL), RLUP,RLDN
  
***    WORK SPACE. 
  
      LOGICAL  MONOTON
      REAL     X,FSC(NSL),DFDLNS(0:NSL1),DLNS(NSL),PRES(LA,NSL) 
      INTEGER  I,K,L,N,INTVL, TOPS,BOTS,INCS, TOPP,BOTP,INCP

      LOGICAL           INFO
      COMMON  /ZZVERBO/ INFO

      EXTERNAL MONVERT,XIT

*---------------------------------------------------------------------
***    CHECK THAT PRLOG IS MONOTONIC.

      CALL MONVERT( PRLOG,NPL, TOPP,BOTP,INCP, MONOTON )

      IF (.NOT.MONOTON)                                        THEN
          IF (INFO) WRITE(6,6000) 'PRESSURE'
          CALL                                     XIT('  Gapl  ',-1 )
      END IF

***    LOOP OVER ALL HORIZONTAL POINTS.
  
      DO  I=1,LA 
  
***        COMPUTE LOCAL SIGMA VALUES AND INVERSE OF DELTA LN SIGMA. 
  
          DO  L=1,NSL
              FSC(L)    = PRES(I,L)-PSLOG(I)
              PRES(I,L) = FSC(L)
          END DO
  
          IF (I.EQ.1)                                          THEN

***            CHECK THAT THE LOCAL SIGMA VECTOR IS MONOTONIC
***            (AND IN WHAT DIRECTION). DO THIS ONLY ONCE.

              CALL MONVERT( FSC,NSL, TOPS,BOTS,INCS, MONOTON )

              IF (.NOT.MONOTON)                                THEN
                  IF (INFO) WRITE(6,6000) 'LOCAL SIGMA'
                  CALL                             XIT('  Gapl  ',-2 )
              END IF

          END IF

          DO  L=TOPS,BOTS-INCS,INCS
              DLNS(L) = 1. / (PRES(I,L+INCS)-PRES(I,L))
          END DO
  
***        GET A COLUMN OF FP (ON PRESSURE LEVELS).
  
          DO  L=TOPS,BOTS,INCS
              FSC(L) = FS(I,L) 
          END DO
  
***        COMPUTE VERTICAL DERIVATIVE OVER ALL PRESSURE INTERVALS.
  
          DO  L=TOPS,BOTS-INCS,INCS
              DFDLNS(L+INCS) = (FSC(L+INCS)-FSC(L)) *DLNS(L)
          END DO
          DFDLNS(TOPS     ) = RLUP
          DFDLNS(BOTS+INCS) = RLDN
  
***        LOOP OVER PRESSURE LEVELS TO BE INTERPOLATED. X IS THE 
***        LN(SIGMA)=LN(PRES)-LN(PS) VALUE OF REQUIRED PRESSURE LEVEL.
  
          K = TOPS
          DO  N=TOPP,BOTP,INCP

              X = PRLOG(N)-PSLOG(I) 
  
***            FIND WHICH SIGMA INTERVAL WE ARE IN.
  
              DO  L=K,BOTS,INCS
                  INTVL = L 
                  IF (X.LT.PRES(I,L)) GOTO 100 
              END DO

              INTVL = BOTS +INCS
  100         K     = INTVL-INCS

              IF (K.EQ.TOPS-INCS) K = TOPS
  
***            NOW INTERPOLATE AT THIS POINT.
  
              FP(I,N)  = FSC(K)+DFDLNS(INTVL)*(X-PRES(I,K))

          END DO
  
      END DO
  
      RETURN
*-----------------------------------------------------------------------

 6000 FORMAT(' Gapl: ',A,' levels are not monotonous.'/)

      END 
      SUBROUTINE INTPZ (FRP,FRXP,FRF,FRX,NIS,NI,RLUP,RLDN)

      IMPLICIT none

      INTEGER  NIS,NI
      REAL     RLUP,RLDN
      REAL     FRF(NIS),FRX(NIS)
      REAL     FRP,FRXP

***    JAN 25/94 - M.GIGUERE. (BASEE SUR LA ROUTINE VERTINT2 DU MC2)
***    MODIFICATION POUR INSERER RLUP ET RLDN 

***    THIS PROGRAM DETERMINES THE VALUE FRP AT POINT FRXP BY CUBIC
***    INTERPOLATION OF THE FUNCTION FRF. THIS FUNCTION IS KNOWN AT NI
***    POINTS AND TO EACH POINT CORRESPONDS THE COORDINATE FRX. THE POINTS
***    MUST BE STORED IN SUCH A WAY THAT THE COORDINATES FRX ARE IN A
***    MONOTONICALLY INCREASING ORDER.

***    THIS SUBROUTINE HANDLES ANY POINT LOCATED OUTSIDE THE GRID OR
***    LOCATED BETWEEN THE BOUNDARY AND THE FIRST INTERIOR POINT.

***    AUTHOR    ANDRE ROBERT                             FEB  1980

      INTEGER  I,PNIA,PNNIM
      REAL     PRDA,PRDB,PRXD,PRSAF,PRSBF,PRSAD,PRSBD

*-----------------------------------------------------------
      PNNIM = NI-1

      IF (FRXP.LE.FRX(1))                                      THEN

***        ON EST PLUS BAS QUE LA PREMIERE COUCHE,
***        EXTRAPOLE AVEC RLDN.

          FRP = FRF(1)-RLDN*(FRX(1)-FRXP)

          GOTO 99

      ELSE IF (FRXP.GE.FRX(NI))                                THEN

***        ON EST AU-DESSUS DE LA PLUS HAUTE COUCHE
***        EXTRAPOLE AVEC RLUP.

          FRP = FRF(NI)+RLUP*(FRXP-FRX(NI))

          GOTO 99

      ELSE IF (FRXP.LE.FRX(2))                                 THEN

***        INTERPOLATION LINEAIRE DANS LA PREMIERE COUCHE.

          PRXD = (FRXP-FRX(1))/(FRX(2)-FRX(1))
          FRP  = (1.0-PRXD)*FRF(1)+PRXD*FRF(2)

          GOTO 99

      ELSE IF (FRXP.GE.FRX(PNNIM))                             THEN

***        INTERPOLATION LINEAIRE DANS LA DERNIERE COUCHE.

          PRXD = (FRXP-FRX(PNNIM))/(FRX(NI)-FRX(PNNIM))
          FRP  = (1.0-PRXD)*FRF(PNNIM)+PRXD*FRF(NI)

          GOTO 99

      ELSE

***        INTERPOLATION CUBIQUE ENTRE LA BASE DE LA PLUS HAUTE COUCHE
***        ET LE SOMMET DE LA PLUS BASSE COUCHE.

***        DETERMINE LA BASE.

          DO  I=2,PNNIM
              IF (FRXP.GE.FRX(I)) PNIA = I
          END DO

***        THE NEXT STEP CONSISTS IN PREPARING THE INFORMATION REQUIRED
***        BY THE INTERPOLATION SUBROUTINE. THE DERIVATIVE IS
***        CALCULATED AT POINTS A AND B
***        IT MUST BE NOTED THAT IT REQUIRES A REALLESS VALUE
***        OF FRXP.

          PRXD = (FRXP-FRX(PNIA))/(FRX(PNIA+1)-FRX(PNIA))
          PRDA = ((FRF(PNIA+1)-FRF(PNIA-1))/(FRX(PNIA+1)-FRX(PNIA-1)))*
     +           (FRX(PNIA+1)-FRX(PNIA))
          PRDB = ((FRF(PNIA+2)-FRF(PNIA))/(FRX(PNIA+2)-FRX(PNIA)))*
     +           (FRX(PNIA+1)-FRX(PNIA))

***        FITS A CUBIC TO THE VALUES FRF(PNIA) AND FRF(PNIA+1) OF THE
***        FUNCTION FRF AT POINTS PNIA AND PNIA+1. THIS CUBIC ALSO FITS
***        THE DERIVATIVES PRDA AND PRDB OF THE FUNCTION AT BOTH POINTS.
***        THE VALUE OF FRP IS CALCULATED AT POINT FRX AND RETURNED TO
***        THE CALLING PROGRAM.

***        THE VALUE FRX MUST VARY FROM 0 TO 1 FROM POINT A TO
***        POINT B. THE FOUR CUBIC SPLINES PRSAF,PRSBF,PRSAD AND PRSBD
***        ARE USED FOR THE INTERPOLATION.

***        ERREUR SI LE TEST SUIVANT EST VRAI.

          IF (PRXD.LT.0.0 .OR. PRXD.GT.1.0) 
     +        CALL                                 XIT('  Intpz ',-1 )

          PRSAF = (1.0+2.0*PRXD)*(1.0-PRXD)*(1.0-PRXD)
          PRSBF = (3.0-2.0*PRXD)*     PRXD *     PRXD
          PRSAD = (1.0    -PRXD)*(1.0-PRXD)*     PRXD
          PRSBD = (1.0    -PRXD)*     PRXD *     PRXD

***        RESULTAT

          FRP = FRF(PNIA)*PRSAF+FRF(PNIA+1)*PRSBF+PRDA*PRSAD-PRDB*PRSBD

      END IF

 99   CONTINUE

      RETURN
*---------------------------------------------------------------------

      END
      SUBROUTINE PLAGC (GALCFLD,PLFLD,PHI,ZGALC,WK1,WK2,IJ,MAXLEV,
     +                  MAXPSIG,NWDS,NSL,NPSL,RLUP,RLDN,GRAV)

      IMPLICIT none

***    JAN 25/94 - M.GIGUERE.

***    INTERPOLATION, EN Z,  DU CHAMP PLFLD AUX NIVEAUX PHI/GRAV
***    A GALCFLD AUX NIVEAUX ZGALC.

      INTEGER  NWDS,NPSL,IJ,MAXPSIG,NSL,MAXLEV
      REAL     GALCFLD(NWDS,NPSL),ZGALC(IJ,MAXPSIG)
      REAL     PLFLD(NWDS,NSL),PHI(NWDS,NSL)
      REAL     WK1(NSL),WK2(NSL)
      REAL     RLUP,RLDN,GRAV

      INTEGER  I,K,LVL

*-----------------------------------------------------------------------
      DO  I=1,NWDS

***        INVERTION DES COLONNES DE PLFLD ET PHI CAR LES VALEURS
***        DE PHI DOIVENT ETRE MONOTONES CROISSANTES.

          DO  K=1,NSL
              WK1(K) = PLFLD(I,NSL-K+1)
              WK2(K) = PHI(I,NSL-K+1)/GRAV
          END DO


          DO  K=1,NPSL
              LVL = NPSL-K+1
              CALL INTPZ( GALCFLD(I,LVL),ZGALC(I,LVL),
     +                    WK1,WK2,NSL,NSL,RLUP,RLDN )
          END DO

***        POUR L'ENSEMBLE DES COLONNES

      END DO

      RETURN
*-----------------------------------------------------------------------

      END
      SUBROUTINE SCALFAC (SF,F,NI,NJ,XP,YP,RS,SCF,COR)

      IMPLICIT none

***    JUL 07/93 - S.BINER, S.TURNER, B.DENIS

***    CETTE SUBROUTINE FAIT LE CALCUL DU FACTEUR D'ECHELLE, SF,
***    ET DU FACTEUR DE CORIOLIS, F.

***    NOTE: LE FACTEUR D'ECHELLE EST CELUI QU'ON SYMBOLISE
***          HABITUELLEMENT PAR M, ET NON PAS PAR S.

      INTEGER  NI,NJ
      REAL     SF(NI,NJ),F(NI,NJ),XP,YP,RS,SCF,COR

      INTEGER  I,J
      REAL     X,Y,SN

*-----------------------------------------------------------------
      DO  J=1,INT(NJ)
          Y = YP-FLOAT( J )
          DO  I=1,INT(NI)
              X       = XP-FLOAT( I )
              SF(I,J) = X**2 + Y**2
              SN      = (RS  - SF(I,J))/(RS + SF(I,J))

              SF(I,J) = SCF  / (1.0 + SN)
              F(I,J)  = COR  * SN
          END DO
      END DO

      RETURN
*--------------------------------------------------------------------------

      END
      SUBROUTINE SDET (COSD,SIND,MNTH,DAY)

      IMPLICIT none

***    APR 10/81 - J.D.HENDERSON 

*     THIS ROUTINE CALCULATES THE SOLAR DECLINATION FOR A GIVEN DAY.
*     ------------------------------------------------------------------

*     D  I  C  T  I  O  N  A  R  Y     (REFER TO SDET AND SDET1). 

*       COSD   - COSINE OF SOLAR DECLINATION. 
*       DAY    - DAY COUNTER -1. 
*       DAYPYR - DAYS IN YEAR.
*       DEC    - 23.5*PI/180*COS(2*PI*(DY-173)/365),SOLAR DECLINATION. 
*       DECMAX - 23.5*PI/180, MAXIMUM SOLAR DECLINATION(RADIANS). 
*       DY     - DAY COUNTER. 
*       EQNX   - EQUINOX,22 JUNE (=173).
*       JDYACC - VARIABLE FOR DAY OF MONTH DETERMINATION. 
*       JY     - LATITUDINAL GRID-POINT INDEX.
*       LVAL   - MONTH INDEX. 
*       MAXDAY - MAXIMUM ALLOWED DAY IN YEAR. 
*       MNTHDY - IDENTIFICATION FOR DAY OF MONTH. 
*       MONTH  - DAYS IN EACH MONTH, BEGINNING WITH JANUARY. 
*       MNTH   - MONTH OF THE YEAR FOR GIVEN DAY.(JAN=1, FEB=2,...,DEC=12)
*       SEASON - (DY-173)/365, TIME PARAMETER IN SOLAR DECLINATION. 
*       SIND   - SINE OF SOLAR DECLINATION. 

*     NOTE THAT THE THE YEAR IS ALWAYS 365 DAYS LONG.

*     ------------------------------------------------------------------
      REAL     COSD,DAY,DAYPYR,DEC,DECMAX,DY,EQNX,SEASON,SIND
      INTEGER  JDY,JDYACC,JYR,L,LVAL,MAXDAY,MNTH,MNTHDY,MONTH(12)
      REAL*8   PI

*     ------------------------------------------------------------------
      DATA     MONTH / 31,28,31,30,31,30,31,31,30,31,30,31 /
*     ------------------------------------------------------------------

      DATA     DAYPYR / 365. /,
     +         DECMAX / 0.4101523743 /,
     +         EQNX   / 173. / 

*-------------------------------------------------------------------- 
      PI     = 4.0*ATAN( 1.D0 )

      MAXDAY = DAYPYR + 1.0E-2
      JDY    = DAY    + 1.0E-2
      JYR    = JDY/MAXDAY
      JDY    = JDY-JYR*MAXDAY
      JDYACC = 0

      DO  L=1,12 
          LVAL   = L
          JDYACC = JDYACC+MONTH(L)
          IF (JDY.LE.JDYACC) GOTO 100 
      END DO
  100 MNTHDY = MONTH(LVAL)-JDYACC+JDY 

      MNTH   = LVAL 
      DY     = JDY
      SEASON = (DY-EQNX)/DAYPYR 

*     ------------------------------------------------------------------
*     EQNX      = JUNE 22
*     APIHELION = JULY 1
*     ------------------------------------------------------------------

      DEC    = DECMAX*COS( 2.0*PI*SEASON ) 
      SIND   = SIN( DEC ) 
      COSD   = COS( DEC ) 

      RETURN
*---------------------------------------------------------------------

      END 
      SUBROUTINE XFILT (FF,F,NI,NJ,NK,IS,JS,KS,CF)

      IMPLICIT none

***    ANDRE ROBERT.     SEPTEMBER 1981.

***    THIS PROGRAM APPLIES A GENERAL SYMMETRIC THREE POINT
***    OPERATOR TO F ALONG THE X-AXIS AND STORES THE RESULT
***    IN FF. IF DESIRED, THE RESULT CAN BE STORED IN THE
***    ORIGINAL FIELD F. VALUES ARE GENERATED AT ALL POINTS.

***    FF=F+CF*FXX

      INTEGER  NI,NJ,NK,IS,JS,KS
      REAL     FF(NI,NJ,NK),F(NI,NJ,NK),CF

      INTEGER  IT,ITT,JT,KT, I,J,K
      REAL     TEMP,FXX

*---------------------------------------------------------------------
      IT  = NI-IS-1
      ITT = IT+1
      JT  = NJ-JS
      KT  = NK-KS

      DO  J=1,JT
          DO  K=1,KT
              TEMP           = F(1,J,K)
              FF(1,J,K)      = TEMP
              DO  I=2,IT
                  FXX       = F(I+1,J,K)+TEMP-2.0*F(I,J,K)
                  TEMP      = F(I,J,K)
                  FF(I,J,K) = F(I,J,K)+CF*FXX
              END DO
              FF(ITT,J,K)   = F(ITT,J,K)
          END DO
      END DO

      RETURN
*---------------------------------------------------------------------

      END
      SUBROUTINE YFILT (FF,F,NI,NJ,NK,IS,JS,KS,CF)

      IMPLICIT none

***    ANDRE ROBERT.     SEPTEMBER 1981.

***    THIS PROGRAM APPLIES A GENERAL SYMMETRIC THREE POINT
***    OPERATOR TO F ALONG THE Y-AXIS AND STORES THE RESULT
***    IN FF. IF DESIRED, THE RESULT CAN BE STORED IN THE
***    ORIGINAL FIELD F. VALUES ARE GENERATED AT ALL POINTS.

***    FF=F+CF*FYY

      INTEGER  NI,NJ,NK,IS,JS,KS
      REAL     FF(NI,NJ,NK),F(NI,NJ,NK),CF

      INTEGER  IT,JT,JTT,KT, I,J,K
      REAL     TEMP,FYY

*---------------------------------------------------------------------
      IT  = NI-IS
      JT  = NJ-JS-1
      JTT = JT+1
      KT  = NK-KS

      DO  I=1,IT
          DO  K=1,KT
              TEMP          = F(I,1,K)
              FF(I,1,K)     = TEMP
              DO  J=2,JT
                  FYY       = F(I,J+1,K)+TEMP-2.0*F(I,J,K)
                  TEMP      = F(I,J,K)
                  FF(I,J,K) = F(I,J,K)+CF*FYY
              END DO
              FF(I,JTT,K)   = F(I,JTT,K)
          END DO
      END DO

      RETURN
*---------------------------------------------------------------------

      END
