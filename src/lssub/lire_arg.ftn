#     if !defined (nombre_de_parametres)
#         define   nombre_de_parametres 91
#     endif
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
      subroutine lire_arg ( GCLE,GDEF,GNAM,NBRGEN,
     +                      IONAM,IPOS,NBRUNT )

      IMPLICIT        none

***    Reads NBRGEN general and any number of local io parametres
***    as well as a maximum of NBRUNT+1 positional arguments.
***    The actual number of such arguments is returned in IPOS.

***    AUTHOR : AVP 20/07 - B.Dugas.

***    Revision 2019/12/05 13:13  dugas
***     Faire appel a DIAG_CCARD plutot qu'a CCARD

***    Input/Output arguments.

      INTEGER         NBRGEN,IPOS,NBRUNT
      CHARACTER*512   GDEF(NBRGEN),IONAM(NBRUNT+1),
     +                GNAM(NBRGEN)
      CHARACTER*16    GCLE(NBRGEN)

***    Local io parametres.

      INTEGER         NBRPAR
      PARAMETER     ( NBRPAR = nombre_de_parametres )
      CHARACTER*512   LDEF(NBRPAR),
     +                LNAM(NBRPAR)
      CHARACTER*16    LCLE(NBRPAR)

      DATA lcle(01)/ 'LON'    /,lnam(01)/ '****' /,ldef(01)/ '  -1' /,
     +     lcle(02)/ 'LON'    /,lnam(02)/ '****' /,ldef(02)/ '****' /,
     +     lcle(03)/ 'LAT'    /,lnam(03)/ '****' /,ldef(03)/ '  -1' /,
     +     lcle(04)/ 'LAT'    /,lnam(04)/ '****' /,ldef(04)/ '****' /,
     +     lcle(05)/ 'LRT'    /,lnam(05)/ '****' /,ldef(05)/ '****' /,
     +     lcle(06)/ 'LMT'    /,lnam(06)/ '****' /,ldef(06)/ '****' /,
     +     lcle(07)/ 'KTR'    /,lnam(07)/ '****' /,ldef(07)/ '****' /,
     +     lcle(08)/ 'KUV'    /,lnam(08)/ '****' /,ldef(08)/ '   1' /,
     +     lcle(09)/ 'LRLMT'  /,lnam(09)/ '****' /,ldef(09)/ '****' /,
     +     lcle(10)/ 'NPG'    /,lnam(10)/ '****' /,ldef(10)/ '****' /,
     +     lcle(11)/ 'NPG'    /,lnam(11)/ '****' /,ldef(11)/ '****' /

      DATA lcle(12)/ 'T1'     /,lnam(12)/ '****' /,ldef(12)/ '  -1' /,
     +     lcle(13)/ 'T1'     /,lnam(13)/ '****' /,ldef(13)/ '****' /,
     +     lcle(14)/ 'T1'     /,lnam(14)/ '****' /,ldef(14)/ '****' /,
     +     lcle(15)/ 'T1'     /,lnam(15)/ '****' /,ldef(15)/ '****' /,
     +     lcle(16)/ 'T1'     /,lnam(16)/ '****' /,ldef(16)/ '****' /,
     +     lcle(17)/ 'T1'     /,lnam(17)/ '****' /,ldef(17)/ '****' /,
     +     lcle(18)/ 'T2'     /,lnam(18)/ '****' /,ldef(18)/ '****' /,
     +     lcle(19)/ 'T2'     /,lnam(19)/ '****' /,ldef(19)/ '****' /,
     +     lcle(20)/ 'T2'     /,lnam(20)/ '****' /,ldef(20)/ '****' /,
     +     lcle(21)/ 'T2'     /,lnam(21)/ '****' /,ldef(21)/ '****' /,
     +     lcle(22)/ 'T2'     /,lnam(22)/ '****' /,ldef(22)/ '****' /

      DATA lcle(23)/ 'T2'     /,lnam(23)/ '****' /,ldef(23)/ '****' /,
     +     lcle(24)/ 'T3'     /,lnam(24)/ '****' /,ldef(24)/ '****' /,
     +     lcle(25)/ 'T3'     /,lnam(25)/ '****' /,ldef(25)/ '****' /,
     +     lcle(26)/ 'T3'     /,lnam(26)/ '****' /,ldef(26)/ '****' /,
     +     lcle(27)/ 'T3'     /,lnam(27)/ '****' /,ldef(27)/ '****' /,
     +     lcle(28)/ 'T3'     /,lnam(28)/ '****' /,ldef(28)/ '****' /,
     +     lcle(29)/ 'T3'     /,lnam(29)/ '****' /,ldef(29)/ '****' /,
     +     lcle(30)/ 'DELT'   /,lnam(30)/ '****' /,ldef(30)/ '****' /,
     +     lcle(31)/ 'A'      /,lnam(31)/ '****' /,ldef(31)/ '****' /,
     +     lcle(32)/ 'B'      /,lnam(32)/ '****' /,ldef(32)/ '****' /,
     +     lcle(33)/ 'C'      /,lnam(33)/ '****' /,ldef(33)/ '****' /

      DATA lcle(34)/ 'NAME'   /,lnam(34)/ '****' /,ldef(34)/ '****' /,
     +     lcle(35)/ 'NAME'   /,lnam(35)/ '****' /,ldef(35)/ '****' /,
     +     lcle(36)/ 'NAME'   /,lnam(36)/ '****' /,ldef(36)/ '****' /,
     +     lcle(37)/ 'NAME'   /,lnam(37)/ '****' /,ldef(37)/ '****' /,
     +     lcle(38)/ 'NAME'   /,lnam(38)/ '****' /,ldef(38)/ '****' /,
     +     lcle(39)/ 'NAME'   /,lnam(39)/ '****' /,ldef(39)/ '****' /,
     +     lcle(40)/ 'NAME'   /,lnam(40)/ '****' /,ldef(40)/ '****' /,
     +     lcle(41)/ 'NAME'   /,lnam(41)/ '****' /,ldef(41)/ '****' /,
     +     lcle(42)/ 'NAME'   /,lnam(42)/ '****' /,ldef(42)/ '****' /,
     +     lcle(43)/ 'NAME'    /,lnam(43)/ '****' /,ldef(43)/ '****' /,
     +     lcle(44)/ 'PLV'     /,lnam(44)/ '****' /,ldef(44)/ '****' /


      DATA lcle(45)/ 'PLV'    /,lnam(45)/ '****' /,ldef(45)/ '****' /,
     +     lcle(46)/ 'PLV'    /,lnam(46)/ '****' /,ldef(46)/ '****' /,
     +     lcle(47)/ 'LV1'    /,lnam(47)/ '****' /,ldef(47)/ '  -1' /,
     +     lcle(48)/ 'LV1'    /,lnam(48)/ '****' /,ldef(48)/ '****' /,
     +     lcle(49)/ 'LV1'    /,lnam(49)/ '****' /,ldef(49)/ '****' /,
     +     lcle(50)/ 'LV2'    /,lnam(50)/ '****' /,ldef(50)/ '****' /,
     +     lcle(51)/ 'DLAT1'  /,lnam(51)/ '****' /,ldef(51)/ '****' /,
     +     lcle(52)/ 'DLAT2'  /,lnam(52)/ '****' /,ldef(52)/ '****' /,
     +     lcle(53)/ 'DLON1'  /,lnam(53)/ '****' /,ldef(53)/ '****' /,
     +     lcle(54)/ 'DLON2'  /,lnam(54)/ '****' /,ldef(54)/ '****' /,
     +     lcle(55)/ 'DGRW'   /,lnam(55)/ '****' /,ldef(55)/ '****' /

      DATA lcle(56)/ 'LX'     /,lnam(56)/ '****' /,ldef(56)/ '****' /,
     +     lcle(57)/ 'LY'     /,lnam(57)/ '****' /,ldef(57)/ '****' /,
     +     lcle(58)/ 'NHEM'   /,lnam(58)/ '****' /,ldef(58)/ '****' /,
     +     lcle(59)/ 'NHEM'   /,lnam(59)/ '****' /,ldef(59)/ '****' /,
     +     lcle(60)/ 'NINTYP' /,lnam(60)/ '****' /,ldef(60)/ '****' /,
     +     lcle(61)/ 'NOUTYP' /,lnam(61)/ '****' /,ldef(61)/ '****' /,
     +     lcle(62)/ 'DEF'    /,lnam(62)/ '****' /,ldef(62)/ '  -1' /,
     +     lcle(63)/ 'D'      /,lnam(63)/ '****' /,ldef(63)/ '****' /,
     +     lcle(64)/ 'SCAL'   /,lnam(64)/ '****' /,ldef(64)/ '****' /,
     +     lcle(65)/ 'I'      /,lnam(65)/ '****' /,ldef(65)/ '****' /,
     +     lcle(66)/ 'I'      /,lnam(66)/ '****' /,ldef(66)/ '****' /

      DATA lcle(67)/ 'I'      /,lnam(67)/ '****' /,ldef(67)/ '****' /,
     +     lcle(68)/ 'J'      /,lnam(68)/ '****' /,ldef(68)/ '****' /,
     +     lcle(69)/ 'J'      /,lnam(69)/ '****' /,ldef(69)/ '****' /,
     +     lcle(70)/ 'J'      /,lnam(70)/ '****' /,ldef(70)/ '****' /,
     +     lcle(71)/ 'K'      /,lnam(71)/ '****' /,ldef(71)/ '****' /,
     +     lcle(72)/ 'K'      /,lnam(72)/ '****' /,ldef(72)/ '****' /,
     +     lcle(73)/ 'K'      /,lnam(73)/ '****' /,ldef(73)/ '****' /,
     +     lcle(74)/ 'L'      /,lnam(74)/ '****' /,ldef(74)/ '****' /,
     +     lcle(75)/ 'L'      /,lnam(75)/ '****' /,ldef(75)/ '****' /,
     +     lcle(76)/ 'L'      /,lnam(76)/ '****' /,ldef(76)/ '****' /,
     +     lcle(77)/ 'M'      /,lnam(77)/ '****' /,ldef(77)/ '****' /

      DATA lcle(78)/ 'M'      /,lnam(78)/ '****' /,ldef(78)/ '****' /,
     +     lcle(79)/ 'M'      /,lnam(79)/ '****' /,ldef(79)/ '****' /,
     +     lcle(80)/ 'N'      /,lnam(80)/ '****' /,ldef(80)/ '****' /,
     +     lcle(81)/ 'N'      /,lnam(81)/ '****' /,ldef(81)/ '****' /,
     +     lcle(82)/ 'N'      /,lnam(82)/ '****' /,ldef(82)/ '****' /,
     +     lcle(83)/ 'LABEL'  /,lnam(83)/ '****' /,ldef(83)/ '****' /,
     +     lcle(84)/ 'LABEL'  /,lnam(84)/ '****' /,ldef(84)/ '****' /,
     +     lcle(85)/ 'LABEL'  /,lnam(85)/ '****' /,ldef(85)/ '****' /,
     +     lcle(86)/ 'D60'    /,lnam(86)/ '****' /,ldef(86)/ '****' /,
     +     lcle(87)/ 'KIND'   /,lnam(87)/ '****' /,ldef(87)/ '****' /,
     +     lcle(88)/ 'KIND'   /,lnam(88)/ '****' /,ldef(88)/ '****' /

      DATA lcle(89)/ 'KIND'   /,lnam(89)/ '****' /,ldef(89)/ '****' /,
     +     lcle(90)/ 'KIND'   /,lnam(90)/ '****' /,ldef(90)/ '****' /,
     +     lcle(91)/ 'KIND'   /,lnam(91)/ '****' /,ldef(91)/ '****' /

***    Local work fields.

      INTEGER         I,NBRCLE
      CHARACTER*4     DEF_PKTYP
      CHARACTER*16,   DIMENSION(:), ALLOCATABLE :: CLES
      CHARACTER*512,  DIMENSION(:), ALLOCATABLE :: NAM,DEF
      CHARACTER*512   PBLOC(NBRPAR)
      LOGICAL         PVAL

      COMMON         /PARBLOC/ PVAL,PBLOC
      COMMON         /ZZDEFPK/ DEF_PKTYP

*----------------------------------------------------------------------
      DEF_PKTYP = 'PK84'

***    Allocate CCARD work fields.

      NBRCLE = NBRGEN + NBRPAR + NBRUNT + 1
      ALLOCATE ( CLES(NBRCLE),NAM(NBRCLE),DEF(NBRCLE) )

***    DEFINE LOGICAL SWITCHES FOUND IN COMMON.

      PVAL = .false.

***    Copy the NBRGEN general parametres.

      DO  I=1,NBRGEN
          CLES(I) = GCLE(I)
          NAM (I) = GNAM(I)
          DEF (I) = GDEF(I)
      END DO

***    Copy the NBRPAR local parametres.

      DO  I=1,NBRPAR
          CLES(NBRGEN+I) = LCLE(I)
          NAM (NBRGEN+I) = LNAM(I)
          DEF (NBRGEN+I) = LDEF(I)
      END DO

***    Define positional ccard character parameters.

      DO  I=NBRGEN+NBRPAR+1,NBRCLE
          CLES(I) = '-.'
          NAM (I) = '  '
          DEF (I) = '  '
      END DO

***    RETRIEVE PARAMETERS.

      CALL DIAG_CCARD( CLES,DEF,NAM,NBRCLE, IPOS )

***    Save the NBRGEN general parametres.

      DO  I=1,NBRGEN
          GNAM(I) = NAM(I)
      END DO

***    Check for non-generic parametre
***    bloc values while filling PBLOC.

      DO  I=1,NBRPAR
          PBLOC(I) = NAM(NBRGEN+I)
          IF (NAM(NBRGEN+I).NE.'****') PVAL = .true.
      END DO

***    Save the IPOS positional parametres.

      DO  I=1,MIN( IPOS,NBRUNT+1 )
          IONAM(I) = NAM(NBRGEN+NBRPAR+I)
      END DO

      DEALLOCATE ( CLES,NAM,DEF )

      RETURN
*----------------------------------------------------------------------

      END
      LOGICAL FUNCTION rpbloc (KEYWORD,VALUE)

***    Auteur: B.Dugas, RPN - 25 novembre 1992.

***    Cette fonction consulte l'etat du bloc de parametre pouvant
***    etre decode par JCLPNT et peut retourner la valeur d'un
***    parametre s'il a ete defini a l'appel du programme.

      INTEGER       NBRPAR
      PARAMETER   ( NBRPAR = nombre_de_parametres )

      CHARACTER*(*) KEYWORD,VALUE

      CHARACTER*512 PBLOC(NBRPAR)
      LOGICAL       PVAL

      COMMON       /PARBLOC/ PVAL,PBLOC

      INTEGER       NBR
      CHARACTER*16  NOMLOC_RPBLOC(NBRPAR)
      SAVE          NOMLOC_RPBLOC

      DATA  NOMLOC_RPBLOC     /
     +     'LON'    , 'LON2'  , 'LAT'    , 'LAT2'   , 'LRT'    ,
     +     'LMT'    , 'KTR'   , 'KUV'    , 'LRLMT'  , 'NPG'    ,
     +     'NPG2'   , 'T1'    , 'T12'    , 'T13'    , 'T14'    ,
     +     'T15'    , 'T16'   , 'T2'     , 'T22'    , 'T23'    ,
     +     'T24'    , 'T25'   , 'T26'    , 'T3'     , 'T32'    ,
     +     'T33'    , 'T34'   , 'T35'    , 'T36'    , 'DELT'   ,
     +     'A'      , 'B'     , 'C'      , 'NAME'   , 'NAME2'  ,
     +     'NAME3'  , 'NAME4' , 'NAME5'  , 'NAME6'  , 'NAME7'  ,
     +     'NAME8'  , 'NAME9' , 'NAME10' , 'PLV'    , 'PLV2'   ,
     +     'PLV3'   , 'LV1'   , 'LV12'   , 'LV13'   , 'LV2'    ,
     +     'DLAT1'  , 'DLAT2' , 'DLON1'  , 'DLON2'  , 'DGRW'   ,
     +     'LX'     , 'LY'    , 'NHEM'   , 'NHEM2'  , 'NINTYP' ,
     +     'NOUTYP' , 'DEFAUT', 'D'      , 'SCAL'   , 'I'      ,
     +     'I2'     , 'I3'    , 'J'      , 'J2'     , 'J3'     ,

     +     'K'      , 'K2'    , 'K3'     , 'L'      , 'L2'     ,
     +     'L3'     , 'M'     , 'M2'     , 'M3'     , 'N'      ,
     +     'N2'     , 'N3'    , 'LABEL'  , 'LABEL2' , 'LABEL3' ,
     +     'D60'    , 'KIND'  , 'KIND2'  , 'KIND3'  , 'KIND4'  ,
     +     'KIND5'            /

*----------------------------------------------------------------------
      rpbloc = .FALSE.
      VALUE  = ' '

      IF (KEYWORD.EQ.' ')                                      THEN

          rpbloc = PVAL

      ELSE IF (PVAL)                                           THEN

          DO  NBR = 1 , NBRPAR

              IF (KEYWORD.NE.NOMLOC_RPBLOC(NBR)) CYCLE

              IF (PBLOC(NBR).NE.'S/O'  .AND.
     +            PBLOC(NBR).NE.'N/A'  .AND.
     +            PBLOC(NBR).NE.'****')                        THEN

                  VALUE  = PBLOC(NBR)
                  rpbloc = .TRUE.

              END IF

              EXIT

          END DO

          IF (NBR.EQ.NBRPAR+1) VALUE = '****'

      END IF
      
      RETURN
*----------------------------------------------------------------------

      END
