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
C     $Log: lowio92.ftn,v $
C     Revision 3.10  2014/12/03 23:26:35  dugas
C     Enlever les enonces EQUIVALENCE.
C
C     Revision 3.9  2014/10/16 12:00:43  dugas
C     Modifications commandees par la compilation avec GFORTRAN et IFORT.
C
C     Revision 3.8  2014/09/25 18:42:03  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.7  2008/04/25 20:51:40  dugas
C     Corriger l'usage des macros pour r.gppf.
C
C     Revision 3.6  2003/10/24 21:05:48  dugas
C     Implementer du code compatible RS6000
C
C     Revision 3.5  2000/03/17 03:17:48  armnrbd
C     Isoler le code CRAY et NEC64 par des macros specifiques
C     dans les routines BF1BI64, BF1I64B, CF1CI64 et CF1I64C.
C
C     Revision 3.4  1999/06/23 20:52:18  armnrbd
C     Remplacer un appel a DABORT par un appel a XIT dans DECODR2.
C
C     Revision 3.3  1999/04/08 19:55:13  armnrbd
C     Utiliser le comdeck MACHTYPE.CDK.
C     Tenir compte de BIGENDI dans ENCODR2, DECODR2 et GETMACH.
C
C     Revision 3.2  1999/01/19 20:09:08  armnrbd
C     Renommer ABORT en DABORT (conflit avec la version systeme).
C
C     Revision 3.1  1994/11/17  14:13:41  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:55:48  13:55:48  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.1  94/06/21  23:19:07  armnrbd
C     Ajouter la routine GETMACH.
C     
C     Revision 2.0  93/10/13  13:31:55  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.2  93/03/08  13:12:03  13:12:03  armnrbd (Bernard Dugas)
C     Corriger une erreur de position de declaration.
C     
C     Revision 1.1  92/11/12  14:08:00  armnrbd
C     Troisieme iteration PKTYP.
C     
C     Revision 1.0  92/11/11  16:02:43  armnrbd
C     Initial revision
C     
      SUBROUTINE GETMACH( IBM,CRAY,IEEE )

***    ESSENTIALLY BASED ON CODE BY J. STACEY - APR 15/92

***    INTEGER REPRESENTATIONS FOR THE HIGH ORDER 32-BITS OF SQRT(2)
***    ARE DIFFERENT ON IBM, CRAY AND IEEE MACHINES.  WE USE THIS TO
***    TEST FOR THE EXACT TYPE OF MACHINE.

      IMPLICIT      none

      LOGICAL       IBM,CRAY,IEEE

      INTEGER       IAMIBM,IAMCRAY,IAMIEEE
      PARAMETER  (  IAMIBM =1092001950 )
      PARAMETER  (  IAMCRAY=1073853700 )
      PARAMETER  (  IAMIEEE=1073127582 )

      REAL*8        A,B
      INTEGER       I,J,K
      POINTER     ( PI,I(1) ),( PJ,J(1) ),( PK,K(2) )

#     include      "machtype.cdk"

      REAL          DATUM(2)
      INTEGER      IDATUM(2)
      EQUIVALENCE ( DATUM,IDATUM )
*------------------------------------------------------------------------------
!!!   EQUIVALENCE ( A,I ),( B,J ),( A,K )
      PI = LOC( A ) ; PJ = LOC( B ) ; PK = PI

***    INITIALIZATION.

      IBM  = .FALSE.
      CRAY = .FALSE.
      IEEE = .FALSE.

***    ISOLATE THE MOST SIGNIFICANT 32-BITS OF SQRT(2).

      A    = SQRT( 2.0D0 )

      IF (BIGENDI.EQ.1 .OR.
     +   (BIGENDI.EQ.0 .AND. MACHINE*INTSIZE.EQ.1))            THEN

          J(1)    = I(1)
          IF (A .EQ. B) J(1) = ISHFT(J(1),-32)

      ELSE

          J(1) = K(2)

      END IF

***    TEST FOR IBM, CRAY OR IEEE MACHINE.

      IBM  = J(1) .EQ. IAMIBM
      CRAY = J(1) .EQ. IAMCRAY
      IEEE = J(1) .EQ. IAMIEEE

      RETURN
      END
      SUBROUTINE BF1BI64(BF,IEEEF,NF)

***    JUL 21/92 - J. STACEY - CORRECT SHIFT OF LEAST SIGNIFICANT BITS.
***    MAY 22/92 - J. STACEY - REMOVE ASSUMPTION THAT INTEGER SIZE IS
***                            NECESSARILY THE SAME AS THE REAL SIZE.
***    APR 15/92 - J. STACEY -

***    CONVERT FLOATING POINT: IBM 64-BIT TO IEEE 64-BIT.

      IMPLICIT     INTEGER (A-Z)

      INTEGER      BF(*),IEEEF(*)
      REAL*8,      SAVE :: ZERO=0.0
      INTEGER      IZERO
      POINTER     (PIZ,IZERO(1))
      INTEGER      IEXP(15),IMAN(15)           !IEEE NORMALIZATION FACTORS.

#     include     "machtype.cdk"

      DATA         IEXP/-4,-3,-3,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1/
      DATA         IMAN/ 0,-1,-1,-2,-2,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3/
*------------------------------------------------------------------------
      PIZ = LOC( ZERO ) !!!   EQUIVALENCE (ZERO,IZERO)

      MASK2 = IBITS(-1,0,7)               !7 BIT MASK FOR IBM EXPONENT.
      MASK4 = IBITS(-1,0,4)               !4 BIT MASK FOR FIRST HEX DIGIT.
      MASK5 = IBITS(-1,0,11)              !11 BIT MASK FOR IEEE EXPONENT.
#     if defined (NEC64)
      IF (INTSIZE .EQ. 1) THEN            !64' INTEGERS/64' REALS.
      MASK1 = IBSET(0,63)                 !1 BIT MASK FOR SIGN BIT.
      MASK3 = IBITS(-1,0,52)              !52 BIT MASK FOR IEEE MANTISSA.
      DO 100 I=1,NF
        IF (BF(I) .NE. IZERO(1)) THEN
            J = IAND(ISHFT(BF(I),-52),MASK4)
            IEEEF(I)  = IOR(IOR(
     1        IAND(BF(I),MASK1),
     2        ISHFT(IAND(4*(IAND(ISHFT(BF(I),-56),MASK2)-64)
     3              +1023+IEXP(J),MASK5),52)
     4                         ),
     5        IAND(ISHFT(BF(I),IMAN(J)),MASK3)
     6                     )
        ELSE
            IEEEF(I) = 0
        END IF
  100 CONTINUE
      END IF
#     else
      IF (INTSIZE .EQ. 2) THEN            !32' INTEGERS/64' REALS.
      MASK1 = IBSET(0,31)                 !1 BIT MASK FOR SIGN BIT.
      MASK3 = IBITS(-1,0,20)              !20 BIT MASK FOR IEEE MANTISSA.
      DO 101 I=1,NF*2,2
        IF (BF(I) .NE. IZERO(1)) THEN
            J = IAND(ISHFT(BF(I),-20),MASK4)
            IEEEF(I)  = IOR(IOR(
     1        IAND(BF(I),MASK1),
     2        ISHFT(IAND(4*(IAND(ISHFT(BF(I),-24),MASK2)-64)
     3              +1023+IEXP(J),MASK5),20)
     4                         ),
     5        IAND(ISHFT(BF(I),IMAN(J)),MASK3)
     6                     )
CBUG        IEEEF(I+1) = ISHFT(BF(I+1),IMAN(J))
            IEEEF(I+1) = IOR(
     1                       ISHFT(BF(I),32+IMAN(J)),
     2                       ISHFT(BF(I+1),IMAN(J))
     3                      ) 
        ELSE
            IEEEF(I) = 0
            IEEEF(I+1) = 0
        END IF
  101 CONTINUE
      END IF
#     endif
      RETURN
      END
      SUBROUTINE BF1I64B(IEEEF,BF,NF)

***    JUL 23/92 - J. STACEY - CORRECT NORMALIZATION OF IBM NUMBERS.
***    JUL 21/92 - J. STACEY - CORRECT SHIFT OF LEAST SIGNIFICANT BITS.
***    MAY 22/92 - J. STACEY - REMOVE ASSUMPTION THAT INTEGER ARRAYS
***                            ARE 64-BITS IN NATIVE FORMAT.
***    APR 15/92 - J. STACEY -
***   
***    CONVERT FLOATING POINT, IEEE 64-BIT TO IBM 64-BIT.  THIS ROUTINE
***    ASSUMES THAT BOTH "BF" AND "IEEEF" ARE 64-BITS IN NATIVE MODE.

      IMPLICIT     INTEGER (A-Z)

      INTEGER      BF(*),IEEEF(*)
      INTEGER      IZERO
      REAL*8,      SAVE :: ZERO=0.0
      LOGICAL      ERROR
      POINTER     (PIZ,IZERO(1))

#     include     "machtype.cdk"

      LOGICAL          INFO
      COMMON /ZZVERBO/ INFO
*--------------------------------------------------------------------
      PIZ = LOC( ZERO ) !!!   EQUIVALENCE (ZERO,IZERO)

      MASK2 = IBITS(-1,0,7)               !7 BIT MASK FOR IBM EXPONENT.
      MASK4 = IBITS(-1,0,4)               !4 BIT MASK FOR FIRST HEX DIGIT.
      MASK5 = IBITS(-1,0,11)              !11 BIT MASK FOR IEEE EXPONENT.
#     if defined (NEC64)
      IF (INTSIZE .EQ. 1) THEN            !64-BIT MACHINE.
      MASK1 = IBSET(0,63)                 !1 BIT MASK FOR SIGN BIT.
      MASK3 = IBITS(-1,0,52)              !52 BIT MASK FOR IEEE MANTISSA.
      IMPLIED = IBSET(0,52)               !IMPLIED ON BIT IN IEEE MANTISSA.
      DO 100 I=1,NF
        IF (IEEEF(I) .NE. 0) THEN
          IBMEXP = IAND(ISHFT(IEEEF(I),-52),MASK5) - 1023 + 1
          IMOVE  = MOD(IBMEXP,4)
          IF (IMOVE .LE. 0) THEN
          IMOVE  = IMOVE + 4
          IBMEXP = IBMEXP/4 + 64
          ELSE
          IBMEXP = IBMEXP/4 + 64 + 1
          ENDIF
          BF(I)  = IOR(IOR(
     1             IAND(IEEEF(I),MASK1),
     2             ISHFT(IAND(IBMEXP,MASK2),56)
     3                     ),
     4             ISHFT(IAND(IEEEF(I),MASK3)+IMPLIED,IMOVE-1)
     5                     )
        ELSE
            BF(I)  = IZERO(1)
        END IF
  100 CONTINUE
      END IF
#     else
      IF (INTSIZE .EQ. 2) THEN            !32-BIT INTEGERS WITH 64-BIT REALS.
      MASK1 = IBSET(0,31)                 !1 BIT MASK FOR SIGN BIT.
      MASK3 = IBITS(-1,0,20)              !20 BIT MASK FOR IEEE MANTISSA.
      IMPLIED = IBSET(0,20)               !IMPLIED ON BIT IN IEEE MANTISSA.
      DO 101 I=1,NF*2,2
        IF (IEEEF(I) .NE. 0 .OR. IEEEF(I+1) .NE. 0) THEN
          IBMEXP = IAND(ISHFT(IEEEF(I),-20),MASK5) - 1023 + 1
          IMOVE  = MOD(IBMEXP,4)
	  IF (IMOVE .LE. 0) THEN
	  IMOVE = IMOVE + 4
	  IBMEXP = IBMEXP/4 + 64
	  ELSE
          IBMEXP = IBMEXP/4 + 64 + 1
	  ENDIF
	  MASK0  = IBITS(-1,0,IMOVE)
          ITEMP  = IEEEF(I+1)
CBUG      BF(I)  = IOR(IOR(
CBUG 1             IAND(IEEEF(I),MASK1),
CBUG 2             ISHFT(IAND(IBMEXP,MASK2),24)
CBUG 3                     ),
CBUG 4             ISHFT(IAND(IEEEF(I),MASK3)+IMPLIED,IMOVE)
CBUG 5                 )
          BF(I)  = IOR(IOR(IOR(
     1             IAND(IEEEF(I),MASK1),
     2             ISHFT(IAND(IBMEXP,MASK2),24)
     3                         ),
     4             ISHFT(IAND(IEEEF(I),MASK3)+IMPLIED,IMOVE-1)
     5                     ),
     6             IAND(ISHFT(ITEMP,-32+IMOVE-1),MASK0)
     7                 )
          BF(I+1) = IAND(ISHFT(ITEMP,IMOVE-1),ISHFT(-1,IMOVE-1))
        ELSE
          BF(I)   = IZERO(1)
          BF(I+1) = 0
        END IF
  101 CONTINUE
      END IF
#     endif
***    HANDLE 0.0, UNDERFLOW AND OVERFLOW :
***   
***    IF THE EXPONENT FIELD GOES NEGATIVE THEN THE IEEE NUMBER WAS
***    EITHER 0.0 OR TOO SMALL TO REPRESENT IN IBM, IN EITHER CASE
***    SET THE IBM RESULT TO 0.0.
***   
***    IF THE EXPONENT FIELD OVERFLOWS, CALL DABORT.

      ERROR = .FALSE.
      IF (INTSIZE .EQ. 1) THEN
        JSHFT = -52
      ELSE
        JSHFT = -20
      END IF
      DO 200 I=1,NF*INTSIZE,INTSIZE
          IBMEXP = IAND(ISHFT(IEEEF(I),JSHFT),MASK5) - 1023 + 1
	  IMOVE = MOD ( IBMEXP,4)
          IF (IMOVE .LE. 0) THEN
          IMOVE  = IMOVE + 4
          IBMEXP = IBMEXP/4 + 64
          ELSE
          IBMEXP = IBMEXP/4 + 64 + 1
          ENDIF
          IF (IBMEXP .LT. 0) THEN
            BF(I) = IZERO(1)
            IF (INTSIZE .EQ. 2) BF(I+1) = 0
          ELSE IF (IBMEXP .GT. MASK2) THEN
            ERROR = .TRUE.
          END IF
  200 CONTINUE
 
      IF (.NOT. ERROR) RETURN
      IF (INFO) WRITE(6,6100)
 6100 FORMAT('0 *ERROR* BF1I64B: IEEE EXPONENT OUT OF RANGE FOR IBM.')
      CALL DABORT
      END

      SUBROUTINE BF2BI64(BF,IEEEF,NF)

***    JUL 27/92 - J. STACEY - CORRECTED NORMALIZATION OF IBM NUMBERS.
***    APR 15/92 - J. STACEY -

***    CONVERT FLOATING POINT: IBM 32-BIT TO IEEE 64-BIT.  THIS ROUTINE
***    IMPLICITLY ASSUMES THAT BOTH "BF" AND "IEEEF" ARE 32-BITS IN NATIVE
***    FORMAT.

      IMPLICIT     INTEGER (A-Z)

      INTEGER      BF(*),IEEEF(*)
      REAL(4),     SAVE :: ZERO=0.0
      INTEGER      IZERO
      POINTER     (PIZ,IZERO(1))
      INTEGER      IEXP(15),IMAN(15)           !IEEE NORMALIZATION FACTORS.
      DATA         IEXP/-4,-3,-3,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1/
      DATA         IMAN/ 0,-1,-1,-2,-2,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3/
*------------------------------------------------------------------------
      PIZ = LOC( ZERO ) !!!   EQUIVALENCE (ZERO,IZERO)     
      MASK1        = IBSET(0,31)          !1 BIT MASK FOR SIGN BIT.
      MASK2        = IBITS(-1,0,7)        !7 BIT MASK FOR IBM EXPONENT.
      MASK3        = IBITS(-1,0,20)       !20 BIT MASK FOR IEEE MANTISSA.
      MASK4        = IBITS(-1,0,4)        !4 BIT MASK FOR FIRST HEX DIGIT.
      MASK5        = IBITS(-1,0,11)       !11 BIT MASK FOR IEEE EXPONENT.
      MASK6        = IBITS(-1,0,3)        !3 BIT MASK FOR L.S.B. IEEE MANT.
      K            = 1
      DO 100 I=1,2*NF-1,2
        IF (BF(K) .NE. IZERO(1)) THEN
            J = IAND(ISHFT(BF(K),-20),MASK4)
            IEEEF(I)  = IOR(IOR(
     1        IAND(BF(K),MASK1),
     2        ISHFT(IAND(4*(IAND(ISHFT(BF(K),-24),MASK2)-64)
     3              +1023+IEXP(J),MASK5),20)
     4                     ),
     5        IAND(ISHFT(BF(K),IMAN(J)),MASK3)
     6                 )
            IEEEF(I+1)  = ISHFT(IAND(BF(K),MASK6),29)
        ELSE
            IEEEF(I)    = 0
            IEEEF(I+1)  = 0
        END IF
        K   = K + 1
  100 CONTINUE
      RETURN
      END

      SUBROUTINE BF2I64B(IEEEF,BF,NF)

***    APR 15/92 - J. STACEY -
***   
***    CONVERT FLOATING POINT, IEEE 64-BIT TO IBM 32-BIT.  THIS ROUTINE
***    ASSUMES THAT BOTH "BF" AND "IEEEF" ARE 32-BITS IN NATIVE MODE.

      IMPLICIT      INTEGER (A-Z)

#     include     "machtype.cdk"

      INTEGER       BF(*),IEEEF(*)
      REAL(4),      SAVE :: ZERO=0.0
      INTEGER       IZERO
      POINTER      (PIZ,IZERO(1))
      LOGICAL       ERROR
      LOGICAL          INFO
      COMMON /ZZVERBO/ INFO
*-------------------------------------------------------------------
      PIZ = LOC( ZERO ) !!!   EQUIVALENCE (ZERO,IZERO)
      IF (INTSIZE .EQ. 1) THEN
          IF (INFO) WRITE(6,6200)
          CALL DABORT
      END IF
      JSHFT   = -52
      MASK1   = IBSET(0,31)         !1 BIT MASK FOR SIGN BIT.
      MASK2   = IBITS(-1,0,7)       !7 BIT MASK FOR IBM EXPONENT.
      MASK3   = IBITS(-1,0,20)      !20 BIT MASK FOR IEEE MANTISSA.
      MASK4   = IBITS(-1,0,4)       !4 BIT MASK FOR FIRST HEX DIGIT.
      MASK5   = IBITS(-1,0,11)      !11 BIT MASK FOR IEEE EXPONENT.
      IMPLIED = IBSET(0,20)         !IMPLIED ON BIT IN IEEE MANTISSA.
      J       = 1
      DO 100 I=1,2*NF-1,2
        IF (IEEEF(I).NE.0 .OR. IEEEF(I+1).NE.0) THEN
CBUG        IBMEXP = IAND(ISHFT(IEEEF(I),-20),MASK5) - 1023
CBUG        IMOVE  = MOD(IBMEXP,4)
CBUG        IBMEXP = IBMEXP/4 + 64 + 1
            IBMEXP = IAND(ISHFT(IEEEF(I),JSHFT),MASK5) - 1023 + 1
            IMOVE  = MOD(IBMEXP,4)
            IF (IMOVE .LE. 0) THEN
            IMOVE  = IMOVE + 4
            IBMEXP = IBMEXP/4 + 64
            ELSE
            IBMEXP = IBMEXP/4 + 64 + 1
            ENDIF
            BF(J) = IOR(IOR(IOR(
     1        IAND(IEEEF(I),MASK1),
     2        ISHFT(IAND(IBMEXP,MASK2),24)
     3                         ),
     4        ISHFT(IAND(IEEEF(I),MASK3)+IMPLIED,IMOVE-1)
     5                     ),
     6        ISHFT(IEEEF(I+1),-32+IMOVE-1)
     7                 )
        ELSE
            BF(J) = IZERO(1)
        END IF
        J = J + 1
  100 CONTINUE
 
***    HANDLE 0.0, UNDERFLOW AND OVERFLOW :
***   
***    IF THE EXPONENT FIELD GOES NEGATIVE THEN THE IEEE NUMBER WAS
***    EITHER 0.0 OR TOO SMALL TO REPRESENT IN IBM, IN EITHER CASE
***    SET THE IBM RESULT TO 0.0.
***   
***    IF THE EXPONENT FIELD OVERFLOWS, CALL DABORT.

      ERROR = .FALSE.
      J     = 1
      DO 20 I=1,2*NF-1,2
CBUG      IBMEXP = IAND(ISHFT(IEEEF(I),-20),MASK5) - 1023
CBUG      IBMEXP = IBMEXP/4 + 64 + 1
          IBMEXP = IAND(ISHFT(IEEEF(I),JSHFT),MASK5) - 1023 + 1
          IMOVE  = MOD(IBMEXP,4)
          IF (IMOVE .LE. 0) THEN
            IMOVE  = IMOVE + 4
            IBMEXP = IBMEXP/4 + 64
          ELSE
            IBMEXP = IBMEXP/4 + 64 + 1
          ENDIF
          IF (IBMEXP .LT. 0) THEN
            BF(J) = IZERO(1)
          ELSE IF (IBMEXP .GT. MASK2) THEN
            ERROR = .TRUE.
          END IF
          J = J + 1
20    CONTINUE
 
      IF (.NOT. ERROR) RETURN
      IF (INFO) WRITE(6,6100)
      CALL DABORT

 6100 FORMAT('0 *ERROR* BF2I64B: IEEE EXPONENT OUT OF RANGE FOR IBM.')
 6200 FORMAT('0 *ERROR* BF2I64B: ILLEGAL CALL, NOT IN 32-BIT MODE')

      END

      SUBROUTINE CF1CI64(CF,IEEEF,NF)

***    APR 15/92 - J. STACEY - A SLIGHTLY MODIFIED VERSION OF "CFCI64",
***                            WHICH WILL COMPILE ON BOTH 32 AND 64 BIT
***                            MACHINES (BUT OBVIOUSLY WILL GIVE USABLE
***                            RESULTS ONLY ON CRAYS).
***   
***    CONVERT FLOATING POINT, CRAY TO IEEE 64-BIT
***   
***   
***    INPUT : CF      CRAY FLOATING POINT NUMBERS
***            NF      NUMBER OF ELEMENTS IN CF
***    OUTPUT: IEEEF   IEEE FLOATING POINT NUMBERS
***   
***    FORMAT :
***              SIGN  EXPONENT  MANTISSA
***    IEEE :     1      11        52
***    CRAY :     1      15        48

      IMPLICIT   INTEGER(A-Z)

      INTEGER    IEEEF(*),CF(*)
      LOGICAL    ERROR
      LOGICAL          INFO
      COMMON /ZZVERBO/ INFO
      PARAMETER (CEXPBIAS=16384)      ! CRAY EXPONENT BIAS (40000(OCTAL))
      PARAMETER (IEXPBIAS=1023)       ! IEEE EXPONENT BIAS
*-------------------------------------------------------------------------------
#     if defined (CRAY)
      MASK1 = IBSET(0,63)
      MASK2 = IBITS(-1,0,15)
      MASK3 = IBITS(-1,0,52)
      MASK4 = IBITS(-1,0,11)

***    SET SIGN BIT, EXPONENT AND MANTISSA IN ONE VECTOR LOOP

      DO 10 I=1,NF
          IEEEF(I) = IOR(IOR(
     1  	  (IAND(CF(I),MASK1))
     2                  ,
     3            (ISHFT((IAND(ISHFT(CF(I),-48),MASK2)-
     4               CEXPBIAS+IEXPBIAS-1),52) )             )
     5                  ,
     6            (IAND(ISHFT(CF(I),5),MASK3))
     7                  )
10    CONTINUE

***    HANDLE 0.0, UNDERFLOW AND OVERFLOW :
***   
***    IF THE EXPONENT FIELD GOES NEGATIVE THEN THE CRAY NUMBER WAS
***    EITHER 0.0 OR TOO SMALL TO REPRESENT IN IEEE, IN EITHER CASE
***    SET THE IEEE RESULT TO ALL 0'S WHICH IS THE IEEE REPRESENTATION 
***    FOR 0.0.  IF THE EXPONENT FIELD OVERFLOWS, CALL DABORT.

      ERROR = .FALSE.
      DO 20 I=1,NF
          IF ((IAND(ISHFT(CF(I),-48),MASK2)-CEXPBIAS+IEXPBIAS-1).LT.0)
     1           IEEEF(I)=0
          IF ((IAND(ISHFT(CF(I),-48),MASK2)-CEXPBIAS+IEXPBIAS-1)
     1           .GT.MASK4)  THEN
            ERROR = .TRUE.
          END IF
20    CONTINUE
 
      IF (.NOT. ERROR) RETURN
      IF (INFO) WRITE(6,6100)
 6100 FORMAT('0***ERROR*** CF1CI64:CRAY EXPONENT OUT OF RANGE.')
      CALL DABORT
#     endif
      END

      SUBROUTINE CF1I64C(IEEEF,CF,NF)

***    APR 15/92 - J. STACEY - A SLIGHTLY MODIFIED VERSION OF "CFI64C",
***                            WHICH WILL COMPILE ON BOTH 32 AND 64 BIT
***                            MACHINES (BUT OBVIOUSLY WILL GIVE USABLE
***                            RESULTS ONLY ON CRAYS).
***   
***    CONVERT FLOATING POINT, IEEE 64-BIT, TO CRAY FLOATING POINT
***   
***    INPUT : IEEEF   IEEE FLOATING POINT NUMBERS (DOUBLE PRECISION)
***            NF      NUMBER OF ELEMENTS IN IEEEF
***    OUTPUT: CF      CRAY FLOATING POINT NUMBERS
***    
***    FORMAT :
***              SIGN  EXPONENT  MANTISSA    UNUSED
***    IEEE :     1      11        52
***    CRAY :     1      15        48
 
      IMPLICIT      INTEGER(A-Z)

      INTEGER       IEEEF(*),CF(*)
      PARAMETER   ( CEXPBIAS=16384)      ! CRAY EXPONENT BIAS (40000(OCTAL))
      PARAMETER   ( IEXPBIAS=1023)       ! IEEE EXPONENT BIAS
      REAL(8),      SAVE :: ZERO=0.0
      INTEGER       IZERO
      POINTER      (PIZ,IZERO(1))
*-----------------------------------------------------------------------------
#     if defined (CRAY)
      PIZ = LOC( ZERO ) !!! EQUIVALENCE ( ZERO,IZERO )
      
      MASK1 = IBSET(0,63)
      MASK2 = IBITS(-1,0,11)
      MASK3 = IBITS(-1,0,52)
      IMPLIED=IBSET(0,52)

***    SET SIGN BIT, EXPONENT AND MANTISSA IN ONE VECTOR LOOP

      DO 10 I=1,NF
        IF (IEEEF(I).NE.0) THEN
          CF(I) = IOR(IOR(
     1  	  (IAND(IEEEF(I),MASK1))
     2              ,
     3            (ISHFT((IAND(ISHFT(IEEEF(I),-52),MASK2)-
     4               IEXPBIAS+CEXPBIAS+1),48) )                 )
     5              ,
     6            (ISHFT((IAND(IEEEF(I),MASK3)+
     7               IMPLIED),-5) )
     8              )
        ELSE
          CF(I) = IZERO(1)
        ENDIF
10    CONTINUE
#     endif 
      RETURN
      END

      SUBROUTINE DECODR2(IEEEF,X)

***    APR 15/92 - J. STACEY - REWRITE OF DECODR FUNCTION, CONVERTED TO
***                            SUBROUTINE CALL TO ALLOW READING OF
***                            64-BIT INTEGER CONSTANTS IN "IEEEF" EVEN
***                            ON 32-BIT MACHINES.
***   
***    PURPOSE: THIS ROUTINES PERFORMS THE INVERSE OPERATION TO
***             ENCODR2. A 64-BIT IEEE REAL*8 VALUE STORED IN "IEEEF"
***             IS CONVERTED INTO A NATIVE SINGLE PRECISION FLOATING
***             POINT VALUE.

      IMPLICIT      REAL (A-H,O-Z), INTEGER (I-N)

      INTEGER       IEEEF(*)

***    INTEGER REPRESENTATIONS FOR THE HIGH ORDER 32-BITS OF SQRT(2)
***    ARE DIFFERENT ON IBM, CRAY AND IEEE MACHINES.  WE USE THIS TO
***    TEST FOR THE EXACT TYPE OF MACHINE BEFORE CONVERTING TO THE
***    STANDARD.

      INTEGER       IAMIBM,IAMCRAY,IAMIEEE
      PARAMETER   ( IAMIBM =1092001950 )
      PARAMETER   ( IAMCRAY=1073853700 )
      PARAMETER   ( IAMIEEE=1073127582 )

      REAL*8        A,B
      INTEGER       I,J,K
      POINTER     ( PI,I(1) ),( PJ,J(1) ),( PK,K(2) )

      LOGICAL       IBM,CRAY,IEEE,IM64BIT,INTEST2

      REAL          DATUM(2)
      INTEGER      IDATUM
      POINTER     ( PID,IDATUM(2) )

      LOGICAL       INFO
      COMMON       /ZZVERBO/ INFO

#     include      "machtype.cdk" 

*------------------------------------------------------------------------------
!!!   EQUIVALENCE ( A,I ),( B,J ),( A,K )
      PI = LOC( A ) ; PJ = LOC( B ) ; PK = PI
      PID = LOC( DATUM(1) ) !!! EQUIVALENCE ( DATUM,IDATUM )

***    INITIALIZATION.

      IBM        = .FALSE.
      CRAY       = .FALSE.
      IEEE       = .FALSE.
      IM64BIT    = .FALSE.

***    USE THE STANDARD TRICK TO SEE IF WE'RE A 32 OR 64 BIT MACHINE.

      A          = DSQRT(2.0D0)

      IF (BIGENDI.EQ.1 .OR.
     +   (BIGENDI.EQ.0 .AND. MACHINE*INTSIZE.EQ.1))            THEN

          J(1)    = I(1)
          IF (A .EQ. B) J = ISHFT(J(1),-32)
          IF (A .EQ. B) IM64BIT  = .TRUE.

      ELSE

          J(1) = K(2)

      END IF

***    TEST FOR IBM, CRAY OR IEEE MACHINE.

      IBM       = J(1) .EQ. IAMIBM
      CRAY      = J(1) .EQ. IAMCRAY
      IEEE      = J(1) .EQ. IAMIEEE
      INTEST2   = INTSIZE .EQ. 2

***    NOW BRANCH ON TYPE OF MACHINE.

      IF (IBM  .AND. IM64BIT) THEN
        IDATUM(1) = 0
        CALL BF1I64B(IEEEF,IDATUM,1)
        X       = DATUM(1)
        RETURN
      END IF

      IF (IBM .AND. .NOT.IM64BIT) THEN
        IF (.NOT.INTEST2) THEN !32' INTEGERS/32' REALS.
          IDATUM(1) = 0
          IDATUM(2) = 0
          CALL BF2I64B(IEEEF,IDATUM,1)
          X       = DATUM(1)
        ELSE                   !32' INTEGERS/64' REALS.
          IDATUM(1) = 0
          IDATUM(2) = 0
          CALL BF1I64B(IEEEF,IDATUM,1)
          X       = DATUM(1)
        END IF
        RETURN
      END IF

      IF (CRAY) THEN
        IDATUM(1) = 0
        CALL CF1I64C(IEEEF,IDATUM,1)
        X       = DATUM(1)
        RETURN
      END IF

      IF (IEEE .AND. IM64BIT) THEN
        IDATUM(1)     = IEEEF(1)
        X             = DATUM(1)
        RETURN
      END IF

      IF (IEEE .AND. .NOT.IM64BIT) THEN
        MASK1 = IBSET(0,31)                 !1 BIT MASK FOR SIGN BIT.
        MASK2 = IBITS(-1,0,8)               !8 BIT MASK FOR IEEE 32' EXPONENT.
        MASK3 = IBITS(-1,0,20)              !20 BIT MASK FOR IEEE 32' MANTISSA.
        MASK4 = IBITS(-1,0,3)               !3 BIT MASK FOR 3' OF IEEE MANTISSA.
        MASK5 = IBITS(-1,0,11)              !11 BIT MASK FOR IEEE 64' EXPONENT.
        IF (.NOT.INTEST2) THEN !INTEGERS ARE 32', REALS ARE 32'.
          IF (IEEEF(1).NE.0 .OR. IEEEF(2).NE.0) THEN
            IDATUM(1) = IOR(IOR(IOR(
     1                  IAND(IEEEF(1),MASK1),
     2              ISHFT(IAND(IBITS(IEEEF(1),20,11)-1023+127,MASK2),23)
     3                          ),
     4              ISHFT(IAND(IEEEF(1),MASK3),3)
     5                           ),
     6              IAND(ISHFT(IEEEF(2),-29),MASK4)
     7                    )
            X      = DATUM(1)
          ELSE
            X      = 0.0
          END IF
          IEXP = IBITS(IEEEF(1),20,11)-1023+127
          IF (IEXP .LT. 0) X = 0.0
          IF (IEXP .GT. MASK2) THEN
            IF (INFO) WRITE(6,6050)
 6050       FORMAT('0 *** ERROR *** DECODR2: IEEE 32-BIT OVERFLOW.')
            CALL XIT(' DECODR2',-1 )
          END IF
        ELSE                   !INTEGERS ARE 32', REALS ARE 64'.
          IDATUM(1) = IEEEF(1)
          IDATUM(2) = IEEEF(2)
          X         = DATUM(1)
        END IF
        RETURN
      END IF

***    IF WE GET HERE THEN WE HAVE A FATAL ERROR (WE DON'T KNOW WHAT MACHINE
***    WE'RE ON OR THE INTERNAL REPRESENTATION OF THE MAGIC CONSTANTS IS
***    DIFFERENT.)

      IF (INFO) WRITE(6,6100)
 6100 FORMAT('0 *** ERROR *** DECODR2: UNKNOWN MACHINE FORMAT.')
      CALL DABORT

      END

      SUBROUTINE ENCODR2(IEEEF,X)

***    APR 15/92 - J. STACEY - REWRITE OF ENCODR FUNCTION, CONVERTED TO
***                            SUBROUTINE CALL TO ALLOW WRITING OF 64-
***                            BIT INTEGER CONSTANTS IN "IEEEF" EVEN ON
***                            32-BIT MACHINES.
***   
***    PURPOSE: THE FLOATING POINT CONSTANT PASSED IN "X" IS CONVERTED
***             TO A 64-BIT IEEE REAL*8 VALUE BEFORE BEING PASSED BACK
***             TO THE CALLING ROUTINE EITHER AS 2 CONSECUTIVE 32-BIT
***             INTEGERS IN "IEEEF" (ON 32-BIT MACHINES) OR AS 1 64-BIT
***             INTEGER IN "IEEEF" (ON 64-BIT MACHINES).

      IMPLICIT      REAL (A-H,O-Z), INTEGER (I-N)

      DIMENSION     IEEEF(*)

***    INTEGER REPRESENTATIONS FOR THE HIGH ORDER 32-BITS OF SQRT(2)
***    ARE DIFFERENT ON IBM, CRAY AND IEEE MACHINES.  WE USE THIS TO
***    TEST FOR THE EXACT TYPE OF MACHINE BEFORE CONVERTING TO THE
***    STANDARD.

      INTEGER       IAMIBM,IAMCRAY,IAMIEEE
      PARAMETER   ( IAMIBM =1092001950 )
      PARAMETER   ( IAMCRAY=1073853700 )
      PARAMETER   ( IAMIEEE=1073127582 )

      REAL*8        A,B
      INTEGER       I,J,K
      POINTER     ( PI,I(1) ),( PJ,J(1) ),( PK,K(2) )

      LOGICAL       IBM,CRAY,IEEE,IM64BIT,INTEST2

      DIMENSION     DATUM(2)
      INTEGER      IDATUM
      POINTER     ( PID,IDATUM(2) )

      LOGICAL       INFO
      COMMON       /ZZVERBO/ INFO

#     include      "machtype.cdk" 

*---------------------------------------------------------------------------
!!!   EQUIVALENCE ( A,I ),( B,J ),( A,K )
      PI = LOC( A ) ; PJ = LOC( B ) ; PK = PI
      PID = LOC( DATUM(1) ) !!! EQUIVALENCE ( DATUM,IDATUM )

***    INITIALIZATION.

      IBM        = .FALSE.
      CRAY       = .FALSE.
      IEEE       = .FALSE.
      IM64BIT    = .FALSE.

***    USE THE STANDARD TRICK TO SEE IF WE'RE A 32 OR 64 BIT MACHINE.

      A          = DSQRT(2.0D0)

      IF (BIGENDI.EQ.1 .OR.
     +   (BIGENDI.EQ.0 .AND. MACHINE*INTSIZE.EQ.1))            THEN

          J(1)   = I(1)
          IF (A .EQ. B) J(1) = ISHFT(J(1),-32)
          IF (A .EQ. B) IM64BIT  = .TRUE.

      ELSE

          J(1) = K(2)

      END IF

***    TEST FOR IBM, CRAY OR IEEE MACHINE.

      IBM       = J(1) .EQ. IAMIBM
      CRAY      = J(1) .EQ. IAMCRAY
      IEEE      = J(1) .EQ. IAMIEEE
      INTEST2   = INTSIZE .EQ. 2

***    NOW BRANCH ON TYPE OF MACHINE.

      IF (IBM  .AND. IM64BIT) THEN
        DATUM(1) = X
        IEEEF(1) = 0
        CALL BF1BI64(IDATUM,IEEEF,1)
        RETURN
      END IF

      IF (IBM .AND. .NOT.IM64BIT) THEN
         IF (.NOT.INTEST2) THEN !32' INTEGERS/32' REALS.
            DATUM(1) = X
            IEEEF(1) = 0
            IEEEF(2) = 0
            CALL BF2BI64(IDATUM,IEEEF,1)
         ELSE                   !32' INTEGERS/64' REALS.
            DATUM(1) = X
            IEEEF(1) = 0
            IEEEF(2) = 0
            CALL BF1BI64(IDATUM,IEEEF,1)
         END IF
         RETURN
      END IF

      IF (CRAY) THEN
        DATUM(1) = X
        IEEEF(1) = 0
        CALL CF1CI64(IDATUM,IEEEF,1)
        RETURN
      END IF

      IF (IEEE .AND. IM64BIT) THEN
        DATUM(1) = X
        IEEEF(1)      = IDATUM(1)
        RETURN
      END IF

      IF (IEEE .AND. .NOT.IM64BIT) THEN
        MASK1 = IBSET(0,31)                 !1 BIT MASK FOR SIGN BIT.
        MASK2 = IBITS(-1,0,8)               !8 BIT MASK FOR IEEE 32' EXPONENT.
        MASK3 = IBITS(-1,0,20)              !20 BIT MASK FOR IEEE 32' MANTISSA.
        MASK4 = IBITS(-1,0,3)               !3 BIT MASK FOR 3' OF IEEE MANTISSA.
        MASK5 = IBITS(-1,0,11)              !11 BIT MASK FOR IEEE 64' EXPONENT.
        IF (.NOT.INTEST2) THEN !INTEGERS ARE 32', REALS ARE 32'.
          DATUM(1) = X
          IF (DATUM(1) .NE. 0.0) THEN
            IEEEF(1) = IOR(IOR(
     1              IAND(IDATUM(1),MASK1),
     2              ISHFT(IAND(IBITS(IDATUM(1),23,8)-127+1023,MASK5),20)
     3                        ),
     4              IAND(ISHFT(IDATUM(1),-3),MASK3)
     5                        )
            IEEEF(2) = ISHFT(IAND(IDATUM(1),MASK4),29)
          ELSE
            IEEEF(1)     = 0
            IEEEF(2)     = 0
          END IF
        ELSE                   !INTEGERS ARE 32', REALS ARE 64'.
          DATUM(1)     = X
          IEEEF(1)     = IDATUM(1)
          IEEEF(2)     = IDATUM(2)
        END IF
        RETURN
      END IF

***    IF WE GET HERE THEN WE HAVE A FATAL ERROR (WE DON'T KNOW WHAT MACHINE
***    WE'RE ON).

      IF (INFO) WRITE(6,6100)
 6100 FORMAT('0 *** ERROR *** ENCODR2: UNKNOWN MACHINE FORMAT.')
      CALL DABORT

      END

      SUBROUTINE IEEEPK(NWDS,LCM,IEEEF)

***    APR 15/92 - J. STACEY -
***   
***    PURPOSE: THIS ROUTINE CONVERTS EITHER FROM IBM 32-BIT, IBM 64-BIT
***             CRAY 64-BIT, IEEE 32-BIT OR IEEE 64-BIT DATA INTO IEEE 64-BIT
***             DATA.  

      IMPLICIT      REAL (A-H,O-Z), INTEGER (I-N)

      DIMENSION     LCM(*),IEEEF(*)

***    INTEGER REPRESENTATIONS FOR THE HIGH ORDER 32-BITS OF SQRT(2)
***    ARE DIFFERENT ON IBM, CRAY AND IEEE MACHINES.  WE USE THIS TO TEST
***    FOR THE EXACT TYPE OF MACHINE BEFORE CONVERTING TO THE STANDARD.

      INTEGER       IAMIBM,IAMCRAY,IAMIEEE
      PARAMETER   ( IAMIBM =1092001950 )
      PARAMETER   ( IAMCRAY=1073853700 )
      PARAMETER   ( IAMIEEE=1073127582 )
      REAL*8        A,B
      INTEGER       I0,J0
      POINTER     ( PI,I0(1) ),( PJ,J0(1) )
      LOGICAL       IBM,CRAY,IEEE,IM64BIT,INTEST2
      REAL,         SAVE :: ZERO=0.0
      INTEGER       IZERO
      POINTER     ( PIZ,IZERO(1))

      LOGICAL       INFO
      COMMON       /ZZVERBO/ INFO

#     include      "machtype.cdk" 
*-------------------------------------------------------------------------------

      PI  = LOC( A ) ; PJ = LOC( B ) !!! EQUIVALENCE ( A,I0 ),( B,J0 )
      PIZ = LOC( ZERO ) !!! EQUIVALENCE ( ZERO,IZERO )

***    INITIALIZATION.

      IBM        = .FALSE.
      CRAY       = .FALSE.
      IEEE       = .FALSE.
      IM64BIT    = .FALSE.

***    USE THE STANDARD TRICK TO SEE IF WE'RE A 32 OR 64 BIT MACHINE.

      A          = DSQRT(2.0D0)
      J0(1)      = I0(1)
      IF (A .EQ. B) THEN
        IM64BIT  = .TRUE.
        J0(1)    = ISHFT(J0(1),-32)
      END IF

***    TEST FOR IBM, CRAY OR IEEE MACHINE.

      IBM       = J0(1) .EQ. IAMIBM
      CRAY      = J0(1) .EQ. IAMCRAY
      IEEE      = J0(1) .EQ. IAMIEEE
      INTEST2   = INTSIZE .EQ. 2

***    NOW BRANCH ON TYPE OF MACHINE.

      IF (IBM  .AND. IM64BIT) THEN
        CALL BF1BI64(LCM,IEEEF,NWDS)
        RETURN
      END IF

      IF (IBM .AND. .NOT.IM64BIT) THEN
        IF (.NOT.INTEST2) THEN !32' INTEGERS/32' REALS.
          CALL BF2BI64(LCM,IEEEF,NWDS)
        ELSE                   !32' INTEGERS/64' REALS.
          CALL BF1BI64(LCM,IEEEF,NWDS)
        END IF
        RETURN
      END IF

      IF (CRAY) THEN
        CALL CF1CI64(LCM,IEEEF,NWDS)
        RETURN
      END IF

      IF (IEEE .AND. IM64BIT) THEN
        DO 20 I=1,NWDS
          IEEEF(I) = LCM(I)
   20   CONTINUE
        RETURN
      END IF

      IF (IEEE .AND. .NOT.IM64BIT) THEN
        K = 1
        MASK1 = IBSET(0,31)
        MASK2 = IBITS(-1,0,8)               !8 BIT MASK FOR IEEE 32' EXPONENT.
        MASK3 = IBITS(-1,0,20)              !20 BIT MASK FOR IEEE 32' MANTISSA.
        MASK4 = IBITS(-1,0,3)               !3 BIT MASK FOR 3' OF IEEE MANTISSA.
        MASK5 = IBITS(-1,0,11)              !11 BIT MASK FOR IEEE 64' EXPONENT.
        IF (.NOT.INTEST2) THEN !INTEGERS ARE 32', REALS ARE 32'.
          DO 30 I=1,2*NWDS-1,2
	    IF (LCM(K) .NE. IZERO(1)) THEN
	      IEEEF(I) = IOR(IOR(
     1                   IAND(LCM(K),MASK1),
     2        ISHFT(IAND(IBITS(LCM(K),23,8)-127+1023,MASK5),20)
     3                        ),
     4                   IAND(ISHFT(LCM(K),-3),MASK3)
     5                    )
              IEEEF(I+1) = ISHFT(IAND(LCM(K),MASK4),29)
	    ELSE
              IEEEF(I)     = 0
              IEEEF(I+1)   = 0
            END IF
            K = K + 1
   30     CONTINUE
        ELSE        !INTEGERS ARE 32', REALS ARE 64'.
          DO 40 I=1,2*NWDS
            IEEEF(I) = LCM(I)
   40     CONTINUE
        END IF
        RETURN
      END IF

***    IF WE GET HERE THEN WE HAVE A FATAL ERROR (WE DON'T KNOW WHAT MACHINE
***    WE'RE ON.

      IF (INFO) WRITE(6,6100)
 6100 FORMAT('0 *** ERROR *** IEEEPK: UNKNOWN MACHINE FORMAT.')
      CALL DABORT

      END

      SUBROUTINE IEEEUP(NWDS,LCM,IEEEF)

***    APR 15/92 - J. STACEY -
***   
***    PURPOSE: THIS ROUTINE CONVERTS FROM IEEE 64-BIT REAL*8 TO EITHER 
***             IBM 32-BIT, IBM 64-BIT CRAY 64-BIT, IEEE 32-BIT OR IEEE 64-BIT 
***             DATA.  

      IMPLICIT      REAL (A-H,O-Z), INTEGER (I-N)

      INTEGER       LCM(*),IEEEF(*)

***    INTEGER REPRESENTATIONS FOR THE HIGH ORDER 32-BITS OF SQRT(2)
***    ARE DIFFERENT ON IBM, CRAY AND IEEE MACHINES.  WE USE THIS TO TEST
***    FOR THE EXACT TYPE OF MACHINE BEFORE CONVERTING TO THE STANDARD.

      INTEGER       IAMIBM,IAMCRAY,IAMIEEE
      PARAMETER   ( IAMIBM =1092001950 )
      PARAMETER   ( IAMCRAY=1073853700 )
      PARAMETER   ( IAMIEEE=1073127582 )
      REAL*8        A,B
      INTEGER       I0,J0
      POINTER     ( PI,I0(1) ),( PJ,J0(1) )
      LOGICAL       IBM,CRAY,IEEE,IM64BIT,INTEST2
      REAL,         SAVE :: ZERO=0.0
      INTEGER       IZERO
      POINTER     ( PIZ,IZERO(1))

      LOGICAL       INFO
      COMMON       /ZZVERBO/ INFO

#     include      "machtype.cdk" 
*--------------------------------------------------------------------------

      PI  = LOC( A ) ; PJ = LOC( B ) !!! EQUIVALENCE ( A,I0 ),( B,J0 )
      PIZ = LOC( ZERO ) !!! EQUIVALENCE ( ZERO,IZERO )

***    INITIALIZATION.

      IBM        = .FALSE.
      CRAY       = .FALSE.
      IEEE       = .FALSE.
      IM64BIT    = .FALSE.

***    USE THE STANDARD TRICK TO SEE IF WE'RE A 32 OR 64 BIT MACHINE.

      A          = DSQRT(2.0D0)
      J0(1)      = I0(1)
      IF (A .EQ. B) THEN
        IM64BIT  = .TRUE.
        J0(1)    = ISHFT(J0(1),-32)
      END IF

***    TEST FOR IBM, CRAY OR IEEE MACHINE.

      IBM        = J0(1) .EQ. IAMIBM
      CRAY       = J0(1) .EQ. IAMCRAY
      IEEE       = J0(1) .EQ. IAMIEEE
      INTEST2    = INTSIZE .EQ. 2

***    NOW BRANCH ON TYPE OF MACHINE.

      IF (IBM  .AND. IM64BIT) THEN
        CALL BF1I64B(IEEEF,LCM,NWDS)
        RETURN
      END IF

      IF (IBM .AND. .NOT.IM64BIT) THEN
        IF (.NOT.INTEST2) THEN !32' INTEGERS/32' REALS
          CALL BF2I64B(IEEEF,LCM,NWDS)
        ELSE                   !32' INTEGERS/64' REALS
          CALL BF1I64B(IEEEF,LCM,NWDS)
        END IF
        RETURN
      END IF

      IF (CRAY) THEN
        CALL CF1I64C(IEEEF,LCM,NWDS)
        RETURN
      END IF

      IF (IEEE .AND. IM64BIT) THEN
        DO 20 I=1,NWDS
          LCM(I) = IEEEF(I)
   20   CONTINUE
        RETURN
      END IF

      IF (IEEE .AND. .NOT.IM64BIT) THEN
        MASK1 = IBSET(0,31)                 !1 BIT MASK FOR SIGN BIT.
        MASK2 = IBITS(-1,0,8)               !8 BIT MASK FOR IEEE 32' EXPONENT.
        MASK3 = IBITS(-1,0,20)              !20 BIT MASK FOR IEEE 32' MANTISSA.
        MASK4 = IBITS(-1,0,3)               !3 BIT MASK FOR 3' OF IEEE MANTISSA.
        MASK5 = IBITS(-1,0,11)              !11 BIT MASK FOR IEEE 64' EXPONENT.
        IF (.NOT.INTEST2) THEN !INTEGERS ARE 32', REALS ARE 32'.
          K = 1
          DO 30 I=1,2*NWDS-1,2
            IF (IEEEF(I).NE.0 .OR. IEEEF(I+1).NE.0) THEN
              LCM(K) = IOR(IOR(IOR(
     1                 IAND(IEEEF(I),MASK1),
     2        ISHFT(IAND(IBITS(IEEEF(I),20,11)-1023+127,MASK2),23)
     3                          ),
     4        ISHFT(IAND(IEEEF(I),MASK3),3)
     5                      ),
     6        IAND(ISHFT(IEEEF(I+1),-29),MASK4)
     7                    )
            ELSE
              LCM(K) = IZERO(1)
            END IF
            K = K + 1
   30     CONTINUE
        ELSE        !INTEGERS ARE 32', REALS ARE 64'.
          DO 40 I=1,2*NWDS
            LCM(I) = IEEEF(I)
   40     CONTINUE
        END IF
        RETURN
      END IF

***    IF WE GET HERE THEN WE HAVE A FATAL ERROR (WE DON'T KNOW WHAT MACHINE
***    WE'RE ON.

      IF (INFO) WRITE(6,6100)
 6100 FORMAT('0 *** ERROR *** IEEEUP: UNKNOWN MACHINE FORMAT.')
      CALL DABORT

      END
