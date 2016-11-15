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
C     $Log: dimgt.ftn,v $
C     Revision 3.6  2016/10/26 16:12  dugas
C     - Restaurer les fonctionnalites de CLRLMT qui avaient ete eliminees en 2013.
C       Elles s'averent requises dans le traitement des troncatures T10->T99 avec
C       des fichiers PK84 ou PK92.
C
C     Revision 3.5  2014/10/16 12:00:36  dugas
C     Modifications commandees par la compilation avec GFORTRAN et IFORT.
C
C     Revision 3.4  2014/09/25 18:31:50  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.3  2013/10/08 01:02:55  bernard
C      - Ajouter le support des LRLMT a 9 chiffres (4+4+1)
C      - Simplifier la routine CLRLMT
C
C     Revision 3.2  1997/02/17 03:55:35  armnrbd
C     Modifier le calcul de LRT et LMT dans DIMGT2 pour
C      corriger une erreur lorsque LRLMT=1010x.
C
C     Revision 3.1  1994/11/17  14:13:05  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:55:18  13:55:18  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:31:35  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.2  93/09/16  23:15:24  23:15:24  armnrbd (Bernard Dugas)
C     Declarer/Sauver MESSAGE dans DIMGT.
C     
C     Revision 1.1  93/08/04  14:11:16  14:11:16  armnrbd (Bernard Dugas)
C     Ajouter la routine CLRLMT.
C     Supporter les deux type de lrlmt par le biais de cette routine.
C     
C     Revision 1.0  92/02/21  11:32:04  armnrbd
C     Initial revision

      SUBROUTINE dimgt (LSR,LA,LR,LM,KTR,LRLMT)

***    DEC 10/90 - B.DUGAS, RPN. (LRLMT PASSE DE AABBC A AAABBBC)
***    AUG 13/90 - M.LAZARE. - CALL XIT IF LR.LT.LM. 
***    JUL 16/79 - J.D.HENDERSON 

***    COMPUTES ROW LENGTH INFORMATION FOR SPECTRAL ARRAYS.
***    LSR(1,M),LSR(2,M) = FIRST WORD IN ROW M OF SPECTRAL,ALP ARRAYS. 

***    LRLMT IS A 5, 7 OR 9 DIGIT NUMBER FORMED BY - 100000*LR+10*LM+KTR,
***                                             BY - 10000*LR+10*LM+KTR,
***                                          OR BY - 1000*LR*10*LM+KTR

***     WHERE  LR = LENGTH OF FIRST SPECTRAL ROW.
***            LM = NUMBER OF SPECTRAL ROWS. 
***           KTR = TRUNCATION TYPE. 
***            LA = TOTAL LENGTH OF SPECTRAL ARRAY (COMPLEX WORDS).

***    INVALID TRUNCATION TYPE RETURNS WITH LA=0.

      IMPLICIT    integer (A-Z)

      DIMENSION   LSR(2,*)

      CHARACTER   MESSAGE*10
      SAVE        MESSAGE

      LOGICAL              INFO
      COMMON     /ZZVERBO/ INFO

      INTEGER     CLRLMT
      EXTERNAL    CLRLMT,XIT

*-------------------------------------------------------------------- 
***    DECOMPOSE LRLMT INTO LR,LM,KTR. 

      LR4 = LRLMT/100000
      LM4 = LRLMT/10-10000*LR4
      LR3 = LRLMT/10000
      LM3 = LRLMT/10-1000*LR3 
      LR2 = LRLMT/1000
      LM2 = LRLMT/10-100*LR2

      IF (LR2 > 99 .OR. LM2 > 99)                              THEN
          IF (LR3 > 999 .OR. LM3 > 999)                        THEN
              LR = LR4
              LM = LM4
          ELSE
              LR = LR3
              LM = LM3
          END IF
      ELSE
          LR = LR2
          LM = LM2
      END IF

      IF (LR < LM .OR. LRLMT > 1000000000)                     THEN
          IF (INFO) WRITE(6,6000) LRLMT,LR,LM
          CALL                                     XIT('  Dimgt ',-1 ) 
      END IF

      IF (INFO                                 .AND. 
     +    LRLMT  .NE.CLRLMT(LR,LM,KTR,.TRUE. ) .AND.
     +    LRLMT  .NE.CLRLMT(LR,LM,KTR,.FALSE.) .AND.
     +    MESSAGE.NE.'DEJA DONNE'            )                 THEN
          WRITE(6,6100) LRLMT,CLRLMT(LR,LM,KTR,.TRUE.)
          MESSAGE = 'DEJA DONNE'
      END IF

      KTR = MOD(LRLMT,10) 
      LA  = 0

      IF (KTR.EQ.0) THEN

***        RHOMBOIDAL TRUNCATION (KTR = 0). 

          LSR(1,1) = 1
          LSR(2,1) = 1

          LMP = LM+1
          DO 250 M=2,LMP
              LSR(1,M) = LSR(1,M-1)+LR
              LSR(2,M) = LSR(2,M-1)+LR+1
  250     CONTINUE

          LA = LSR(1,LMP)-1 

      ELSE IF(KTR.EQ.2) THEN

***        TRIANGULAR TRUNCATION (KTR = 2).

          LSR(1,1) = 1
          LSR(2,1) = 1

          LMP = LM+1
          DO 350 M=2,LMP
              LSR(1,M) = LSR(1,M-1)+LR-(M-2)
              LSR(2,M) = LSR(2,M-1)+LR-(M-3)
  350     CONTINUE

          LA = LSR(1,LMP)-1 

      END IF

      RETURN
*-----------------------------------------------------------------------

 6000 FORMAT(' Dimgt reads lrlmt= ',I10,'. Finds lr,lm= ',2I5,'.')
 6100 FORMAT(' LRLMT received/calculated are: ',2I10,'.')

      END 

      SUBROUTINE dimgt2 (LSR,LA,LR,LM,KTR,LRLMT,KIND,CALC)

***    DEC 10/90 - B.DUGAS, RPN. (LRLMT PASSE DE AABBC A AAABBBC)
***    AUG 13/90 - M.LAZARE. - CALL XIT IF LR.LT.LM. 
***    NOV 27/88 - B.DUGAS.

***    COMPUTES ROW LENGTH INFORMATION FOR SPECTRAL ARRAYS FROM VALUES 
***    OF LRLMT AND KIND. THE RESULT ARE RETURNED AS...

***    LSR(1,M),LSR(2,M) = FIRST WORD IN ROW M OF SPECTRAL,ALP ARRAYS. 

***    LRLMT IS A 5, 7 OR 9 DIGIT NUMBER FORMED BY - 100000*LR+10*LM+KTR,
***                                             BY - 10000*LR+10*LM+KTR,
***                                          OR BY - 1000*LR*10*LM+KTR

***    WHERE  LR = LENGTH OF FIRST SPECTRAL ROW. 
***           LM = NUMBER OF SPECTRAL ROWS.
***          KTR = TRUNCATION TYPE.

***    THE FIVE, SEVEN OR NINE DIGIT FORMAT ARE TO BE USED DEPENDING ON
***    THE  RANGE NEEDED TO CORRECTLY CARRY THE SPECTRAL RESOLUTION INFO.
***    (T/R)8 AND LOWER WILL ALWAYS USE THE 5-DIGIT REPRESENTATION, 
***    WHILE (T/R)9 TO (T/R)98 CAN MAKE USE OF EITHER THE 5 OR 7 DIGIT
***    REPRESENTATION. (T/R)99 TO (T/R)998 ARRAYS WILL USE THE 7-DIGIT
***    REPRESENTATION ONLY.

***    KIND DETERMINES THE PROCESSING OPTION IN THE FOLLOWING WAY
***          = 0 ==> GLOBAL COEFFICIENTS INPUT, GLOBAL OUTPUT GRIDS
***          < 0 ==> ANTI-SYMMETRIC COEFFICIENTS, HEMISPHERIC OUTPUT 
***          > 0 ==> SYMMETRIC COEFFICIENTS,        "     "      " 

***    IF CALC.NE.0 LSR IS NOT ACTUALLY WRITTEN BUT ALL THE OTHER 
***    PARAMETERS ARE RETURNED CORRECTLY.

***    LA = TOTAL LENGTH OF SPECTRAL ARRAY (COMPLEX WORDS).

***    NOTE THAT INVALID TRUNCATION  TYPES RETURN  WITH LA = 0. RHOM-
***    BOIDAL AND TRIANGULAR TRUNCATIONS ONLY ARE SUPPORTED. CONSTANT
***    ZERO COEFFICIENTS IN THE KIND.NE.0 CASES ARE SQUEEZED OUT FROM
***    THE LSR(1,*) COUNT. 

      IMPLICIT    integer (A-Z)
  
      DIMENSION   LSR(2,*)

      CHARACTER   MESSAGE*10
      SAVE        MESSAGE

      LOGICAL              INFO
      COMMON     /ZZVERBO/ INFO

      INTEGER     CLRLMT
      EXTERNAL    CLRLMT,XIT

*-------------------------------------------------------------------- 
***    DECOMPOSE LRLMT INTO LR,LM,KTR. 
  
      LR4 = LRLMT/100000
      LM4 = LRLMT/10-10000*LR4
      LR3 = LRLMT/10000
      LM3 = LRLMT/10-1000*LR3 
      LR2 = LRLMT/1000
      LM2 = LRLMT/10-100*LR2

      IF (LR2 > 99 .OR. LM2 > 99)                              THEN
          IF (LR3 > 999 .OR. LM3 > 999)                        THEN
              LR = LR4
              LM = LM4
          ELSE
              LR = LR3
              LM = LM3
          END IF
      ELSE
          LR = LR2
          LM = LM2
      END IF

      IF (LR < LM .OR. LRLMT > 1000000000)                    THEN
          IF (INFO) WRITE(6,6000) LRLMT,LR,LM
          CALL                                     XIT(' Dimgt2 ',-1 ) 
      END IF

      KTR = MOD(LRLMT,10) 
      LMP = LM+1

      IF (INFO                                 .AND. 
     +    LRLMT  .NE.CLRLMT(LR,LM,KTR,.TRUE. ) .AND.
     +    LRLMT  .NE.CLRLMT(LR,LM,KTR,.FALSE.) .AND.
     +    MESSAGE.NE.'DEJA DONNE'            )                 THEN
          WRITE(6,6100) LRLMT,CLRLMT(LR,LM,KTR,.TRUE.)
          MESSAGE = 'DEJA DONNE'
      END IF

      IF (CALC.EQ.0)                                           THEN 
  
          LSR(1,1) = 1
          LSR(2,1) = 1
  
          IF (KTR.EQ.0)                                        THEN 
  
***            RHOMBOIDAL TRUNCATION.
  
              IF (KIND.EQ.0)                                   THEN 
                  DO 100 M=2,LMP
                      LSR(1,M) = LSR(1,M-1)+LR
                      LSR(2,M) = LSR(2,M-1)+LR+1
  100             CONTINUE
              ELSE
                  LSR(1,2) = LR+1 
                  LSR(2,2) = LR+2 
                  IF (KIND.GT.0)                               THEN 
                      LRSYMM = (LR+MOD(LR,2))/2 
                      DO 110 M=3,LMP
                          LSR(1,M) = LSR(1,M-1)+LRSYMM
                          LSR(2,M) = LSR(2,M-1)+LR+1
  110                 CONTINUE
                  ELSE
                      LRUNSY = LR/2 
                      DO 120 M=3,LMP
                          LSR(1,M) = LSR(1,M-1)+LRUNSY
                          LSR(2,M) = LSR(2,M-1)+LR+1
  120                 CONTINUE
                  END IF
              END IF
  
          ELSE IF (KTR.EQ.2)                                   THEN 
  
***            TRIANGULAR TRUNCATION.
  
              IF (KIND.EQ.0)                                   THEN 
                  DO 200 M=2,LMP
                      LSR(1,M) = LSR(1,M-1)+LR-(M-2)
                      LSR(2,M) = LSR(2,M-1)+LR-(M-3)
  200             CONTINUE
              ELSE
                  LSR(1,2) = LR+1 
                  LSR(2,2) = LR+2 
                  IF (KIND.GT.0)                                   THEN 
                      DO 210 M=3,LMP
                          LSR(1,M) = LSR(1,M-1)+(LR-M+3)/2
                          LSR(2,M) = LSR(2,M-1)+LR-(M-3)
  210                 CONTINUE
                  ELSE
                      DO 220 M=3,LMP
                          LSR(1,M) = LSR(1,M-1)+(LR-M+2)/2
                          LSR(2,M) = LSR(2,M-1)+LR-(M-3)
  220                 CONTINUE
                  END IF
              END IF
  
          ELSE
  
              LSR(1,LMP) = 1
  
          END IF
  
          LA = LSR(1,LMP)-1 
  
      ELSE
  
          LSR1 = 1
          LSR2 = 1
  
          IF (KTR.EQ.0)                                        THEN 
  
***            RHOMBOIDAL TRUNCATION.
  
              IF (KIND.EQ.0)                                   THEN 
                  DO 300 M=2,LMP
                      LSR1 = LSR1+LR
                      LSR2 = LSR2+LR+1
  300             CONTINUE
              ELSE
                  LSR1 = LR+1 
                  LSR2 = LR+2 
                  IF (KIND.GT.0)                               THEN 
                      LRSYMM = (LR+MOD(LR,2))/2 
                      DO 310 M=3,LMP
                          LSR1 = LSR1+LRSYMM
                          LSR2 = LSR2+LR+1
  310                 CONTINUE
                  ELSE
                      LRUNSY = LR/2 
                      DO 320 M=3,LMP
                          LSR1 = LSR1+LRUNSY
                          LSR2 = LSR2+LR+1
  320                 CONTINUE
                  END IF
              END IF
  
          ELSE IF (KTR.EQ.2)                                   THEN 
  
***            TRIANGULAR TRUNCATION.
  
              IF (KIND.EQ.0)                                   THEN 
                  DO 400 M=2,LMP
                      LSR1 = LSR1+LR-(M-2)
                      LSR2 = LSR2+LR-(M-3)
  400             CONTINUE
              ELSE
                  LSR1 = LR+1 
                  LSR2 = LR+2 
                  IF (KIND.GT.0)                                   THEN 
                      DO 410 M=3,LMP
                          LSR1 = LSR1+(LR-M+3)/2
                          LSR2 = LSR2+LR-(M-3)
  410                 CONTINUE
                  ELSE
                      DO 420 M=3,LMP
                          LSR1 = LSR1+(LR-M+2)/2
                          LSR2 = LSR2+LR-(M-3)
  420                 CONTINUE
                  END IF
              END IF
  
          END IF
  
          LA = LSR1-1 
  
      END IF
  
      RETURN
*-----------------------------------------------------------------------

 6000 FORMAT(' Dimgt2 reads lrlmt= ',I10,'. Finds lr,lm= ',2I5,'.')
 6100 FORMAT(' LRLMT received/calculated are: ',2I10,'.')

      END 
      INTEGER FUNCTION CLRLMT (LR,LM,KTR,CONDIT)

      IMPLICIT none

***    AUTHOR: B.Dugas, RPN - 17 juillet 1993.

      INTEGER   LR,LM,KTR
      LOGICAL   CONDIT

***    PURPOSE: BUILDS LRLMT FIELD FROM THE LR, LM AND KTR FIELDS. FIVE-,
***             SEVEN- AND NINE-DIGITS LRLMT's ARE SUPPORTED. T10 TO T99
***             USES THE FIVE OR SEVEN DIGITS REPRESENTATION, DEPENDING
***             ON THE VALUE OF CONDIT AND POSSIBLY OF THE 'LRLMT_DIGITS'
***             ENVIRONMENTAL VARIABLE ; T100 TO T999 USES SEVEN DIGITS
***             AND NINE DIGITS ARE USED IN THE T1000 TO T9999 RANGE

      CHARACTER(LEN=256), SAVE :: EVAL=' '
      INTEGER,            SAVE :: DIGITS

      LOGICAL            INFO
      COMMON   /ZZVERBO/ INFO

      EXTERNAL  XIT,GETENVC

*----------------------------------------------------------------------
***    CONSISTENCY CHECKING.

      IF (LR .LT.0 .OR.
     +   (LM .LT.0 .OR. LR .LT.LM) .OR.
     +   (KTR.LT.0 .OR. KTR.GT.9 ) )                           THEN
          IF (INFO) WRITE(6,6001) LR,LM,KTR
          CALL                                     XIT(' Clrlmt',-1 )
      END IF

      IF (LR > 99 .OR. LM > 99)                                THEN
          IF (LR > 999 .OR. LM > 999)                          THEN
              CLRLMT = (MOD( LR,10000 )*10000+MOD( LM,10000 ))*10+KTR
          ELSE
              CLRLMT = (MOD( LR,1000 )*1000+MOD( LM,1000 ))*10+KTR
          END IF
      ELSE IF (LR < 10 .AND. LM < 10)                          THEN
          CLRLMT = (LR*100+LM)*10+KTR
      ELSE
          IF (EVAL.EQ.' ')                                     THEN
              CALL GETENVC( 'LRLMT_DIGITS',EVAL )
              IF (EVAL == ' ') EVAL = ' 7'
              READ(EVAL,'(BN,I2)') DIGITS
              IF (DIGITS /= 5 .AND. DIGITS /= 7)               THEN
                  IF (INFO) WRITE(6,6002) DIGITS
                  CALL                             XIT(' Clrlmt',-2 )
              END IF                 
          END IF
          IF (DIGITS == 5)                                     THEN
              IF (CONDIT)                                      THEN
                  CLRLMT = (LR*100+LM)*10+KTR
              ELSE
                  CLRLMT = (LR*1000+LM)*10+KTR
              END IF
          ELSE
              IF (CONDIT)                                      THEN
                  CLRLMT = (LR*1000+LM)*10+KTR
              ELSE
                  CLRLMT = (LR*100+LM)*10+KTR
              END IF
          END IF
      END IF

      RETURN
*----------------------------------------------------------------------

 6001 FORMAT(' Clrlmt found lr,lm,ktr= ',3I5,'.')
 6002 FORMAT(" Clrlmt found 'LRLMT_DIGITS' = ",I2)

      END
