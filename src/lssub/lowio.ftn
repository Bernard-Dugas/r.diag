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
!     $Log: lowio.ftn,v $
!     Revision 3.10  2015/04/23 21:15:58  dugas
!     Autres modifications pour GFORTRAN.
!
!     Revision 3.9  2015/04/23 20:41:36  dugas
!     Modifier les comparaisons avec ZERO et corriger certains enonces TRANSFER pour AIX.
!
!     Revision 3.8  2014/12/03 23:26:35  dugas
!     Enlever les enonces EQUIVALENCE.
!
!     Revision 3.7  2014/09/25 18:42:03  dugas
!     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
!
!     Revision 3.6  2008/04/25 20:50:45  dugas
!     Remplacer 'defined f77' par defined (F77)' dans ENCODR/DECODR.
!
!     Revision 3.5  2000/04/04 16:48:57  armnrbd
!     Tenir compte de la definition de constantes hexadecimales
!     avec F90 dans les routines ENCODR et DECODR.
!
!     Revision 3.4  1999/06/23 20:53:02  armnrbd
!     Verifier pour la presence de NaN dans DECODR.
!
!     Revision 3.3  1999/04/08 19:53:31  armnrbd
!     Utiliser le comdeck MACHTYPE.CDK.
!     Tenir compte de BIGENDI dans ENCODR et DECODR.
!     Elaborer differamment les boucles principales des
!     routines GBYTES1, GBYTES2, GBYTES3 et GBYTES4.
!
!     Revision 3.2  1997/11/21 20:58:14  armnrbd
!     Assurer l'alignement de Z dans (DE/EN)CODR.
!
!     Revision 3.1  1994/11/17  14:13:40  armnrbd
!     Messages informatifs quand au passage de la version 2.x a 3.1...
!     1) Les espaces en debut des noms de variables de sont plus pertinents.
!     2) Les grilles complexes de type CMPL sont maintenant supportees.
!     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
!     4) Plusieurs nouvelles cles sont disponibles au demarrage.
!
!     Revision 3.0  94/11/17  13:55:46  13:55:46  armnrbd (Bernard Dugas)
!     *** empty log message ***
!     
!     Revision 2.0  93/10/13  13:31:54  armnrbd
!     Premiere version compatible HP-UX.
!     
!     Revision 1.0  92/02/21  11:33:40  11:33:40  armnrbd (Bernard Dugas)
!     Initial revision
!     

      SUBROUTINE encodr (X,Y)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION Y(2)

!     * APRIL 01/90 - B.DUGAS. CREATES AN IEEE-754 VERSION OF THE
!     *                        AES/CCRN PACKER.
!     * NOV 06/90   - B.DUGAS. CONVERT THE INPUT FROM REAL*4 FORMAT.

!     * THIS ROUTINE TAKES  A VALID IEEE-754 FP NUMBER (X) AND
!     * NORMALIZES IT IN SUCH A WAY THAT IT CAN BE THEN BE SENT
!     * FROM ONE COMPUTER TO ANOTHER. THE OUTPUT (Y) CLOSELY
!     * RESEMBLES A CRAY SINGLE PRECISION FP NUMBER.

!     * THE IEEE-754 64-BIT INPUT WORD LOOKS LIKE: 

!     *           +-----------------------------------------------+
!     *           ! SIGN !  BIASED EXPONENT  !     MANTISSA       !
!     *           +-----------------------------------------------+

!     * WHERE ...

!     *       SIGN     = 1 BIT,
!     *       EXPONENT = 11 BITS (BIASED BY 2**10-1) AND 
!     *       MANTISSA = 52 BITS (WITH AN IMPLICIT LEADING 1).

!     * THE NORMALISED 64-BIT OUTPUT NUMBER LOOKS LIKE: 

!     *           +-----------------------------------------------+
!     *           !   BIASED EXPONENT   !    BIASED MANTISSA      !
!     *           +-----------------------------------------------+

!     * WHERE ...

!     *       EXPONENT = 15 BITS (2**14 BIASED EXPONENT) AND
!     *       MANTISSA = 49 BITS (BIASED BY 2**48-1 WHICH IS THE
!     *                           LARGEST POSSIBLE MANTISSA ON A
!     *                           CRAY COMPUTER).

!-------------------------------------------------------------------
      DIMENSION Z(2)
      REAL(8)   Z8
      POINTER (IZ8,Z8(*))

#     include    "machtype.cdk"

#     if defined (F77)
      DATA INME    / X'7FF00000' /,
     +     INBE    / X'78040000' /,
     +     INMM    / X'FFFFF'    /,
     +     OUB1    / X'7FFFFFF'  /,
     +     OUB2    / X'1FFFFF'   /,
     +     OUMM    / X'1F'       /,
     +     ZEROMAN / X'FFFFFFFF' /,
     +     ZEROEXP / X'8000FFFF' /
#     else
      REAL(8), SAVE :: ZERO8 = Z'8000FFFFFFFFFFFF'
      INTEGER, SAVE :: IZ1   = TRANSFER( Z'8000FFFF',1 )
      INTEGER, SAVE :: IZ2   = TRANSFER( Z'FFFFFFFF',1 )
      INTEGER, SAVE :: INME  = TRANSFER( Z'7FF00000',1 )
      INTEGER, SAVE :: INBE  = TRANSFER( Z'78040000',1 )
      INTEGER, SAVE :: INMM  = TRANSFER( Z'000FFFFF',1 ) 
      INTEGER, SAVE :: OUB1  = TRANSFER( Z'07FFFFFF',1 )
      INTEGER, SAVE :: OUB2  = TRANSFER( Z'001FFFFF',1 )
      INTEGER, SAVE :: OUMM  = TRANSFER( Z'0000001F',1 )
#     endif
      INTEGER, SAVE :: ZERO=0, UN=1, TROIS=3, CINQ=5, LAST=31

!     * (TRANSLATING) STATEMENT FUNCTIONS.
#     if defined (HP) || defined (SUN)
      AND( NUM1, NUM2) = IAND   ( NUM1, NUM2)
      OR ( NUM1, NUM2) = IOR    ( NUM1, NUM2)
#     endif
!     ISHL( NUM,  BITS) = ISHFT ( NUM,  BITS)
      ISHC( NUM,  BITS) = ISHFTC( NUM,  BITS, 32)

      EXTERNAL GO4A8

! ------------------------------------------------------------------- 
!     IZ8 = LOC( Y(1) )

      IF (BIGENDI == 1)                                        THEN
          HIW = 1
          LOW = 2
      ELSE
          HIW = 2
          LOW = 1
      END IF

!     * CONVERT TO REAL*8 FORMAT BEFORE CODING.

      CALL GO4A8 (X,Z)

!     * TAKE CARE OF ZERO VALUES.

      IF (       Z(HIW)      == ZERO .AND. 
     +    ISHFT( Z(LOW),UN ) == ZERO      )                    THEN
#         if defined (F77)
          Y(2) = ZEROMAN
          Y(1) = ZEROEXP
#         else
          Y(1) = IZ1 ; Y(2) = IZ2
!         Z8(1) = ZERO8
#         endif
          RETURN
      END IF

!     * ISOLATE EXPONENT E IN RESULT WORD. ACCOUNT FOR IEEE BIAS.
!     * NOTE THAT INBE = 2*(2**14-(2**10-1)+1)).

      E = ISHFT( AND( INME, Z(HIW)) , -TROIS) + INBE

!     * ISOLATE MANTISSA. ISOLATE THE HIGH 27 BITS FROM Z(HIW) AND 
!     * THE  LOW 20 BITS FROM Z(LOW). THEN, ACCOUNTING FOR THE IM-
!     * PLICIT 2**51 BIT GIVES US A 48-BIT MANTISSA.

      M1 = ISHFT( Z(LOW),-CINQ )
      M2 = AND(   Z(HIW), INMM ) + (INMM+UN)

!     * CHECK IF SIGN BIT IS ON. 

      IF (BTEST( Z(HIW),LAST ))                                THEN

!         * NEGATIVE NUMBER CASE. SUBSTRACT THE MANTISSA FROM THE BIAS.

          M1 = OUB1-M1
          M2 = OUB2-M2

      ELSE

!         * POSITIVE NUMBER CASE. ADD THE BIAS TO THE MANTISSA,
!         * TAKING CARE OF A POSSIBLE CARRY BIT FROM M1 TO M2.

          IF (M1 /= ZERO)                                      THEN
              M1 = M1-UN
              M2 = M2+UN+OUB2
          ELSE
              M1 = OUB1
              M2 = OUB2+M2
          END IF
          
      END IF

!     * JOIN THE NORMALIZED EXPONENT AND MANTISSA.

      Y(2) = OR( ISHC( AND( OUMM,M2 ),-CINQ ),M1 )
      Y(1) = OR( ISHFT( M2,-CINQ ),E )

      RETURN

!--------------------------------------------------------------------
      END

      SUBROUTINE decodr (X,Y)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION Y(2)

!     * APRIL 01/90 - B.DUGAS. CREATES AN IEEE-754 VERSION OF THE
!     *                        AES/CCRN PACKER.
!     * NOV 06/90   - B.DUGAS. CONVERT THE OUTPUT TO REAL*4 FORMAT.

!     * THIS SUBPROGRAM DOES THE OPPOSITE WORK OF ENCODR. IT
!     * TAKES A NORMALIZED FP NUMBER (Y) AND RETURNS A IEEE-754 
!     * VALID FORM (X).  NOTE THAT THE INPUT CLOSELY RESEMBLES A
!     * CRAY SINGLE PRECISION FP NUMBER.

!     * THE NORMALIZED 64-BIT INPUT NUMBER LOOKS LIKE: 

!     *           +-----------------------------------------------+
!     *           !   BIASED EXPONENT   !    BIASED MANTISSA      !
!     *           +-----------------------------------------------+

!     * WHERE ...

!     *       EXPONENT = 15 BITS (2**14 BIASED EXPONENT) AND
!     *       MANTISSA = 49 BITS (MANTISSA BIASED BY 2**48-1,
!     *                           WHICH IS THE LARGEST MANTISSA
!     *                           ON A CRAY COMPUTER).

!     * THE 64-BIT OUTPUT F.P. NUMBER WILL LOOK LIKE: 

!     *           +-----------------------------------------------+
!     *           ! SIGN !  BIASED EXPONENT  !     MANTISSA       !
!     *           +-----------------------------------------------+

!     * WHERE ...

!     *       SIGN     = 1 BIT,
!     *       EXPONENT = 11 BITS (BIASED BY 2**10-1) AND
!     *       MANTISSA = 52 BITS (WITH AN IMPLICIT LEADING 1),

!     * I.E. A STANDARD IEEE-754 DOUBLE PRECISION F.P. NUMBER.

!-------------------------------------------------------------------
      DIMENSION Z(2)
      REAL(8)   Z8
      POINTER  (IZ8,Z8(*))

      LOGICAL              INFO
      COMMON     /ZZVERBO/ INFO

#     include    "machtype.cdk"

#     if defined (F77)
      DATA INME    / X'FFFE0000' /,
     +     INBE    / X'78020000' /,
     +     INM1    / X'1FFFF'    /,
     +     INM2    / X'F0000000' /,
     +     OUB1    / X'FFFFFFF'  /,
     +     OUB2    / X'FFFFF'    /,
     +     EINC    / X'100000'   /,
     +     ZEROMAN / X'FFFFFFFF' /,
     +     ZEROEXP / X'8000FFFF' /
#     else
      REAL(8), SAVE :: ZERO8 = Z'8000FFFFFFFFFFFF'
      INTEGER, SAVE :: IZ1   = TRANSFER( Z'8000FFFF',1 )
      INTEGER, SAVE :: IZ2   = TRANSFER( Z'FFFFFFFF',1 )
      INTEGER, SAVE :: INME  = TRANSFER( Z'FFFE0000',1 )
      INTEGER, SAVE :: INBE  = TRANSFER( Z'78020000',1 )
      INTEGER, SAVE :: INM1  = TRANSFER( Z'0001FFFF',1 ) 
      INTEGER, SAVE :: INM2  = TRANSFER( Z'F0000000',1 )
      INTEGER, SAVE :: OUB1  = TRANSFER( Z'0FFFFFFF',1 )
      INTEGER, SAVE :: OUB2  = TRANSFER( Z'000FFFFF',1 )
      INTEGER, SAVE :: EINC  = TRANSFER( Z'00100000',1 )
#     endif
      INTEGER, SAVE :: ZERO=0, UN=1, TROIS=3, QUATRE=4, C28=-28, LAST=31

!     * (TRANSLATING) STATEMENT FUNCTIONS.
#     if defined (HP) || defined (SUN)
      AND( NUM1, NUM2) = IAND   ( NUM1, NUM2)
      OR ( NUM1, NUM2) = IOR    ( NUM1, NUM2)
#     endif
!     ISHL( NUM,  BITS) = ISHFT (NUM,  BITS)

      EXTERNAL    GO8A4,XIT

! ------------------------------------------------------------------- 
!     IZ8 = LOC( Y(1) )

!     * TAKE CARE OF ZERO VALUES.

#     if defined (F77)
      IF (Y(2) == ZEROMAN .AND. Y(1) == ZEROEXP)               THEN
#     else
!     IF (Z8(1) == ZERO8)                                      THEN
      IF (Y(1) == IZ1 .AND. Y(2) == IZ2)                       THEN
#     endif
          X = ZERO
          RETURN
      END IF

!     * ISOLATE EXPONENT E IN RESULT WORD. CORRECT TO IEEE BIAS.
!     * NOTE THAT INBE = 2*(2**14-(2**10-1)).

      E = ISHFT( AND( INME,Y(1) ) - INBE,TROIS )

!     * ISOLATE MANTISSA. ISOLATE THE 32 BITS FROM Y(2) AND THE
!     * LOW 17 BITS FROM Y(1).  THIS GIVES US A 49-BIT MANTISSA.
!     * SEPARATE IN HIGH 21- AND LOW 28-BIT WORDS AND SUBSTRACT 
!     * THE BIAS. IF THE MANTISSA IS SMALLER THAN THE BIAS, THE
!     * SIGN BIT IS TURNED ON.

      M1 = AND( OUB1,Y(2) )

      M3 = AND( INM2,Y(2) )
      M4 = AND( INM1,Y(1) )

      M3 = ISHFT( M3,C28 )
      M4 = ISHFT( M4,QUATRE )

      M2 = OR( M3,M4 )

      IF (M2.LT.OUB2 .OR. (M2 == OUB2 .AND. M1.LT.OUB1))       THEN

!         * NEGATIVE NUMBER CASE. SUBSTRACT THE MANTISSA FROM THE BIAS.

          S  = IBSET( ZERO,LAST )
          M1 = OUB1-M1
          M2 = OUB2-M2

      ELSE

!         * POSITIVE NUMBER CASE.  SUBSTRACT THE BIAS FROM THE MAN-
!         * TISSA, TAKING CARE OF A POSSIBLE CARRY BIT FROM M1 TO M2.

          S = ZERO

          IF (M1.LT.OUB1)                                      THEN
              M1 = M1+UN
              M2 = M2-UN-OUB2
          ELSE
              M1 = ZERO
              M2 = M2-OUB2
          END IF

      END IF

!     * JOIN THE SIGN, EXPONENT AND MANTISSA. CHECK THAT BIT 21
!     * OF M2 IS SET. IF IT IS NOT, SHIFT LEFT UNTIL IT IS.

  100 IF (M2.LE.OUB2)                                          THEN
          E  = E-EINC
          M1 = ISHFT( M1,UN )
          M2 = ISHFT( M2,UN )
          IF (M1.GT.OUB1)                                      THEN
              M1 = M1-OUB1
              M2 = M2+UN
          END IF
          GOTO 100
      END IF

      IF (BIGENDI == 1)                                        THEN
          Z(2) = ISHFT( M1,QUATRE )
          Z(1) = OR( AND( M2,OUB2 ),OR( S,E ) )
      ELSE
          Z(1) = ISHFT( M1,QUATRE )
          Z(2) = OR( AND( M2,OUB2 ),OR( S,E ) )
      END IF

!     * CHECK FOR IEEE-754 64-BITS INFINITE OR NAN VALUES.

      NaN = ISHFT( 1+(AND( 2047,ISHFT( Z(BIGENDI),-20 ) )),-11 )

      IF (NaN == 0)                                            THEN

!         * GO TO SINGLE PRECISION IEEE-754.

          CALL GO8A4 (Z,X)

      ELSE

          IF (INFO) WRITE(6,6000)
          CALL                                     XIT(' Decodr ',-1 )

      END IF

      RETURN

!--------------------------------------------------------------------
 6000 FORMAT(//'0 *** Error *** DECODR: IEEE overflow.'//)

      END 

      SUBROUTINE gbytesb (JPAK,J,BITS,DIM)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION J(DIM), JPAK(*)

!     * APRIL 04/1990 - B.DUGAS (UQAM)

!     * THIS IS A SUBSET OF THE CRAY  SBYTESB ROUTINE (WRITTEN BY 
!     * J.F. FORTIN, CCRN - DECEMBER/1983).  IN THIS VERSION, THE 
!     * OUTPUT HAS TO FILL 32-BIT INTEGERS, IF THE RESULTS ARE TO
!     * REMAIN IDENTICAL WITH THE CRAY VERSION. AS SUCH,  PACKING
!     * DENSITIES OF 2,4,8 AND 16 ARE SUPPORTED. 32 BIT REAL WORD 
!     * ARE NO LONGER CONSIDERED (THEY SHOULD NOT BE PACKED). NON
!     * 32 BIT-FILLING PACKING DENSITIES WORK FINE, THANK YOU.

!     * THIS SUBROUTINE UNPACKS THE MEMORY WORDS ALREADY PACKED 
!     * WITH SBYTESB. THE  PACKED NUMBERS ARE POSITIVE INTEGERS   
!     * CODED ON "BITS" BITS.  THE UNPACKED ARRAY IS J AND WILL 
!     * BE DIMENSIONNED TO "DIM".

!----------------------------------------------------------------------
!     * "BW" IS THE NUMBER OF BITS PER INTEGER WORD.  THIS VALUE 
!     * IS MACHINE  AND IMPLEMENTATION DEPENDENT. FOR EXAMPLE, A
!     * SINGLE PRECISION NUMBER WARRANTS BW = 64 ON A CRAY COMPU-
!     * TER.

      COMMON /MACHSPEC/ BW

      INTEGER, SAVE :: ZERO=0, UN=1, DEUX=2

!     * (TRANSLATING) STATEMENT FUNCTIONS.
#     if defined (HP) || defined (SUN)
      AND( NUM1, NUM2) = IAND   ( NUM1, NUM2)
#     endif
!     ISHL( NUM,  BITS) = ISHFT (NUM,  BITS)

! ---------------------------------------------------------------------
!     * PACKING DENSITY.

      PACK  = BW/BITS

!     * MASK USED TO SELECT WHICH COLUMN OF "BITS" BITS WILL BE
!     * PROCESSED IN JPAK.

      MASQUE = 2**BITS-UN

!     * NUMBER OF BITS TO SHIFT THE POSITIVE INTEGERS INTO THE 
!     * UNPACKED ARRAY J.

      BSKIP  = (PACK-UN)*BITS
      MASQUE = ISHFT( MASQUE, BSKIP)

      DO  POINT = UN,(PACK-UN)

!         * FILL J WITH  ONE COLUMN OF "BITS" BITS OF JPAK BY 
!         * INCREMENT OF "PACK",  AFTER  SHIFTING THE  "BITS"
!         * COLUMNS TRANSFERED BY BSKIP COLUMNS TO THE RIGHT.

          DO  K = POINT,DIM,PACK
              JJ   = K/PACK + 1
              J(K) = ISHFT( AND( MASQUE, JPAK(JJ)), -BSKIP) 
          END DO

          MASQUE   = ISHFT( MASQUE, -BITS)
          BSKIP    = BSKIP-BITS

      END DO


!     * DO THE LAST SET SEPARATELY BECAUSE THERE IS NO SHIFTING.

      DO  K = PACK,DIM,PACK
          JJ   = K/PACK
          J(K) = AND( MASQUE, JPAK(JJ))
      END DO

      RETURN

!--------------------------------------------------------------------
      END 
 
      SUBROUTINE gbytes1 (JPAK,X,BITS,DIM,XSCALI,XMIN)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION JPAK(*)
      REAL(4)   X(*),XSCALI,XMIN

!     * SEPTEMBER 16/1991 - B.DUGAS (RPN)

!     * THIS IS A SUBSET OF THE CRAY  SBYTESB ROUTINE (WRITTEN BY 
!     * J.F. FORTIN, CCRN - DECEMBER/1983).  IN THIS VERSION, THE 
!     * OUTPUT HAS TO FILL 32-BIT INTEGERS, IF THE RESULTS ARE TO
!     * REMAIN  IDENTICAL WITH THE CRAY VERSION.  THIS IS THE 1/1 
!     * PACKING VERSION.

!     * THIS SUBROUTINE UNPACKS THE MEMORY WORDS ALREADY PACKED 
!     * WITH SBYTES1. THE  PACKED NUMBERS ARE POSITIVE INTEGERS   
!     * CODED ON "BITS" BITS.  THE UNPACKED ARRAY IS X AND WILL 
!     * BE DIMENSIONNED TO "DIM".

!----------------------------------------------------------------------
!     * "BW" IS THE NUMBER OF BITS PER INTEGER WORD.  THIS VALUE 
!     * IS MACHINE  AND IMPLEMENTATION DEPENDENT. FOR EXAMPLE, A
!     * SINGLE PRECISION NUMBER WARRANTS BW = 64 ON A CRAY COMPU-
!     * TER.

      COMMON /MACHSPEC/ BW

!     * (TRANSLATING) STATEMENT FUNCTIONS.
#     if defined (HP) || defined (SUN)
      AND( NUM1, NUM2) = IAND   ( NUM1, NUM2)
#     endif
!     ISHL( NUM,  BITS) = ISHFT (NUM,  BITS)

! ---------------------------------------------------------------------
!     * MASK USED TO SELECT WHICH COLUMN OF "BITS" BITS WILL BE
!     * PROCESSED IN JPAK.

      MASQUE = 2**BITS-1

!     * FILL J WITH  ONE COLUMN OF "BITS" BITS OF JPAK BY 
!     * INCREMENT OF "PACK",  AFTER  SHIFTING THE  "BITS"
!     * COLUMNS TRANSFERED BY BSKIP COLUMNS TO THE RIGHT.

      DO  K = 1,DIM

          IX   = AND( MASQUE,JPAK(K) )
          X(K) = IX * XSCALI + XMIN

      END DO

!--------------------------------------------------------------------
      END 
 
      SUBROUTINE gbytes2 (JPAK,X,BITS,DIM,XSCALI,XMIN)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION JPAK(*)
      REAL(4)   X(2,*),XSCALI,XMIN

!     * SEPTEMBER 16/1991 - B.DUGAS (RPN)

!     * THIS IS A SUBSET OF THE CRAY  SBYTESB ROUTINE (WRITTEN BY 
!     * J.F. FORTIN, CCRN - DECEMBER/1983).  IN THIS VERSION, THE 
!     * OUTPUT HAS TO FILL 32-BIT INTEGERS, IF THE RESULTS ARE TO
!     * REMAIN  IDENTICAL WITH THE CRAY VERSION.  THIS IS THE 2/1 
!     * PACKING VERSION.

!     * THIS SUBROUTINE UNPACKS THE MEMORY WORDS ALREADY PACKED 
!     * WITH SBYTES2. THE  PACKED NUMBERS ARE POSITIVE INTEGERS   
!     * CODED ON "BITS" BITS.  THE UNPACKED ARRAY IS X AND WILL 
!     * BE DIMENSIONNED TO "DIM".

!----------------------------------------------------------------------
!     * "BW" IS THE NUMBER OF BITS PER INTEGER WORD.  THIS VALUE 
!     * IS MACHINE  AND IMPLEMENTATION DEPENDENT. FOR EXAMPLE, A
!     * SINGLE PRECISION NUMBER WARRANTS BW = 64 ON A CRAY COMPU-
!     * TER.

      COMMON /MACHSPEC/ BW

!     * (TRANSLATING) STATEMENT FUNCTIONS.
#     if defined (HP) || defined (SUN)
      AND( NUM1, NUM2) = IAND   ( NUM1, NUM2)
#     endif
!     ISHL( NUM,  BITS) = ISHFT (NUM,  BITS)

! ---------------------------------------------------------------------
!     * PACKING DENSITY.

      PACK = 2
      DIM2 = DIM/PACK

!     * MASK USED TO SELECT WHICH COLUMN OF "BITS" BITS WILL BE
!     * PROCESSED IN JPAK.

      MASQUE2 = 2**BITS-1
      MASQUE1 = ISHFT( MASQUE2, BITS )

      BSKIP1  = -BITS

!     * FILL J WITH  ONE COLUMN OF "BITS" BITS OF JPAK BY 
!     * INCREMENT OF "PACK",  AFTER  SHIFTING THE  "BITS"
!     * COLUMNS TRANSFERED BY BSKIP COLUMNS TO THE RIGHT.

      DO  K = 1,DIM2

          IX2    =        AND( MASQUE2,JPAK(K) )
          IX1    = ISHFT( AND( MASQUE1,JPAK(K) ), BSKIP1 )

          X(2,K) = IX2 * XSCALI + XMIN
          X(1,K) = IX1 * XSCALI + XMIN

      END DO

      IF (MOD(DIM,2) /= 0)
     +    X(1,K) = ISHFT(
     +                    AND( MASQUE1,JPAK(K) ), 
     +                    BSKIP1
     +                  )
     +           * XSCALI + XMIN

!--------------------------------------------------------------------
      END 
 
      SUBROUTINE gbytes3 (JPAK,X,BITS,DIM,XSCALI,XMIN)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION JPAK(*)
      REAL(4)   X(3,*),XSCALI,XMIN

!     * SEPTEMBER 16/1991 - B.DUGAS (RPN)

!     * THIS IS A SUBSET OF THE CRAY  SBYTESB ROUTINE (WRITTEN BY 
!     * J.F. FORTIN, CCRN - DECEMBER/1983).  IN THIS VERSION, THE 
!     * OUTPUT HAS TO FILL 32-BIT INTEGERS, IF THE RESULTS ARE TO
!     * REMAIN  IDENTICAL WITH THE CRAY VERSION.  THIS IS THE 3/1 
!     * PACKING VERSION.

!     * THIS SUBROUTINE UNPACKS THE MEMORY WORDS ALREADY PACKED 
!     * WITH SBYTES3. THE  PACKED NUMBERS ARE POSITIVE INTEGERS   
!     * CODED ON "BITS" BITS.  THE UNPACKED ARRAY IS X AND WILL 
!     * BE DIMENSIONNED TO "DIM".

!----------------------------------------------------------------------
!     * "BW" IS THE NUMBER OF BITS PER INTEGER WORD.  THIS VALUE 
!     * IS MACHINE  AND IMPLEMENTATION DEPENDENT. FOR EXAMPLE, A
!     * SINGLE PRECISION NUMBER WARRANTS BW = 64 ON A CRAY COMPU-
!     * TER.

      COMMON /MACHSPEC/ BW

!     * (TRANSLATING) STATEMENT FUNCTIONS.
#     if defined (HP) || defined (SUN)
      AND( NUM1, NUM2) = IAND   ( NUM1, NUM2)
#     endif
!     ISHL( NUM,  BITS) = ISHFT (NUM,  BITS)

! ---------------------------------------------------------------------
!     * PACKING DENSITY.

      PACK = 3
      DIM3 = DIM/PACK

!     * MASK USED TO SELECT WHICH COLUMN OF "BITS" BITS WILL BE
!     * PROCESSED IN JPAK.

      MASQUE3 = 2**BITS-1
      MASQUE2 = ISHFT( MASQUE3, BITS )
      MASQUE1 = ISHFT( MASQUE2, BITS )

      BSKIP2 = -BITS*1
      BSKIP1 = -BITS*2

!     * FILL J WITH  ONE COLUMN OF "BITS" BITS OF JPAK BY 
!     * INCREMENT OF "PACK",  AFTER  SHIFTING THE  "BITS"
!     * COLUMNS TRANSFERED BY BSKIP COLUMNS TO THE RIGHT.

      DO  K = 1,DIM3

          IX3    =        AND( MASQUE3,JPAK(K) )
          IX2    = ISHFT( AND( MASQUE2,JPAK(K) ), BSKIP2 )
          IX1    = ISHFT( AND( MASQUE1,JPAK(K) ), BSKIP1 )

          X(3,K) = IX3 * XSCALI + XMIN
          X(2,K) = IX2 * XSCALI + XMIN
          X(1,K) = IX1 * XSCALI + XMIN

      END DO

      IF (MOD(DIM,3) /= 0)                                     THEN

          IX1    = ISHFT( AND( MASQUE1,JPAK(K) ), BSKIP1 )
          X(1,K) = IX1 * XSCALI + XMIN

          IF (MOD(DIM,3).GT.1)                                 THEN
              IX2    = ISHFT( AND( MASQUE2,JPAK(K) ), BSKIP2 )
              X(2,K) = IX2 * XSCALI + XMIN
          END IF    

      END IF

!--------------------------------------------------------------------
      END 
 
      SUBROUTINE gbytes4 (JPAK,X,BITS,DIM,XSCALI,XMIN)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION JPAK(*)
      REAL(4)   X(4,*),XSCALI,XMIN

!     * SEPTEMBER 16/1991 - B.DUGAS (RPN)

!     * THIS IS A SUBSET OF THE CRAY  SBYTESB ROUTINE (WRITTEN BY 
!     * J.F. FORTIN, CCRN - DECEMBER/1983).  IN THIS VERSION, THE 
!     * OUTPUT HAS TO FILL 32-BIT INTEGERS, IF THE RESULTS ARE TO
!     * REMAIN  IDENTICAL WITH THE CRAY VERSION.  THIS IS THE 4/1 
!     * PACKING VERSION.

!     * THIS SUBROUTINE UNPACKS THE MEMORY WORDS ALREADY PACKED 
!     * WITH SBYTES4. THE  PACKED NUMBERS ARE POSITIVE INTEGERS   
!     * CODED ON "BITS" BITS.  THE UNPACKED ARRAY IS X AND WILL 
!     * BE DIMENSIONNED TO "DIM".

!----------------------------------------------------------------------
!     * "BW" IS THE NUMBER OF BITS PER INTEGER WORD.  THIS VALUE 
!     * IS MACHINE  AND IMPLEMENTATION DEPENDENT. FOR EXAMPLE, A
!     * SINGLE PRECISION NUMBER WARRANTS BW = 64 ON A CRAY COMPU-
!     * TER.

      COMMON /MACHSPEC/ BW

!     * (TRANSLATING) STATEMENT FUNCTIONS.
#     if defined (HP) || defined (SUN)
      AND( NUM1, NUM2) = IAND   ( NUM1, NUM2)
#     endif
!     ISHL( NUM,  BITS) = ISHFT (NUM,  BITS)

! ---------------------------------------------------------------------
!     * PACKING DENSITY.

      PACK = 4
      DIM4 = DIM/PACK

!     * MASK USED TO SELECT WHICH COLUMN OF "BITS" BITS WILL BE
!     * PROCESSED IN JPAK.

      MASQUE4 = 2**BITS-1
      MASQUE3 = ISHFT( MASQUE4, BITS )
      MASQUE2 = ISHFT( MASQUE3, BITS )
      MASQUE1 = ISHFT( MASQUE2, BITS )

      BSKIP3 = -BITS*1
      BSKIP2 = -BITS*2
      BSKIP1 = -BITS*3

!     * FILL J WITH  ONE COLUMN OF "BITS" BITS OF JPAK BY 
!     * INCREMENT OF "PACK",  AFTER  SHIFTING THE  "BITS"
!     * COLUMNS TRANSFERED BY BSKIP COLUMNS TO THE RIGHT.

      DO  K = 1,DIM4

          IX4    =        AND( MASQUE4,JPAK(K) )
          IX3    = ISHFT( AND( MASQUE3,JPAK(K) ), BSKIP3 )
          IX2    = ISHFT( AND( MASQUE2,JPAK(K) ), BSKIP2 )
          IX1    = ISHFT( AND( MASQUE1,JPAK(K) ), BSKIP1 )

          X(4,K) = IX4 * XSCALI + XMIN
          X(3,K) = IX3 * XSCALI + XMIN
          X(2,K) = IX2 * XSCALI + XMIN
          X(1,K) = IX1 * XSCALI + XMIN

      END DO

      IF (MOD(DIM,4) /= 0)                                     THEN

          IX1    = ISHFT( AND( MASQUE1,JPAK(K) ), BSKIP1 )
          X(1,K) = IX1 * XSCALI + XMIN

          IF (MOD(DIM,4).GT.1)                                 THEN

              IX2    = ISHFT( AND( MASQUE2,JPAK(K) ), BSKIP2 )
              X(2,K) = IX2 * XSCALI + XMIN

              IF (MOD(DIM,4).GT.2)                             THEN
                  IX3    = ISHFT( AND( MASQUE3,JPAK(K) ), BSKIP3 )
                  X(3,K) = IX3 * XSCALI + XMIN
              END IF

          END IF

      END IF

!--------------------------------------------------------------------
      END 
 
      SUBROUTINE sbytesb (JPAK,J,BITS,DIM)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION J(DIM),JPAK(*)

!     * APRIL 04/1990 - B.DUGAS (UQAM)

!     * THIS IS A SUBSET OF THE CRAY  SBYTESB ROUTINE (WRITTEN BY 
!     * J.F. FORTIN, CCRN - DECEMBER/1983).  IN THIS VERSION, THE 
!     * OUTPUT HAS TO FILL 32-BIT INTEGERS, IF THE RESULTS ARE TO
!     * REMAIN IDENTICAL WITH THE CRAY VERSION. AS SUCH,  PACKING
!     * DENSITIES OF 2,4,8 AND 16 ARE SUPPORTED. 32 BIT REAL WORD 
!     * ARE NO LONGER CONSIDERED (THEY SHOULD NOT BE PACKED). NON
!     * 32 BIT-FILLING PACKING DENSITIES WORK FINE, THANK YOU.

!     * THIS SUBROUTINE PACKS AN ARRAY J OF  DIMENSION IDIM, CONTAINING
!     * POSITIVE INTEGERS CODED ON NBITS BITS, INTO AN ARRAY JPAK WHICH 
!     * DIMENSION IS <= DIM.

!----------------------------------------------------------------------
!     * "BW" IS THE NUMBER OF BITS PER INTEGER WORD.  THIS VALUE 
!     * IS MACHINE  AND IMPLEMENTATION DEPENDENT. FOR EXAMPLE, A
!     * SINGLE PRECISION NUMBER WARRANTS BW = 64 ON A CRAY COMPU-
!     * TER.

      COMMON /MACHSPEC/ BW

      INTEGER, SAVE :: ZERO=0, UN=1, DEUX=2

!     * (TRANSLATING) STATEMENT FUNCTIONS.
#     if defined (HP) || defined (SUN)
      OR ( NUM1, NUM2) = IOR    ( NUM1, NUM2)
#     endif
!     ISHL( NUM,  BITS) = ISHFT (NUM,  BITS)

! ------------------------------------------------------------------- 
!     * PACKING DENSITY.

      PACK = BW/BITS

!     * REMAINDER OF ELEMENTS WHICH REQUIRE BUT MAY NOT FILL A JPACK
!     * ELEMENT SPACE.

      RMNDR = MOD( DIM, PACK)

!     * DIMENSION OF PACKED ARRAY JPAK. 

      IF (RMNDR == 0)                                          THEN
          PDIM = DIM/PACK
      ELSE
          PDIM = (DIM/PACK)+UN
      END IF

!     * INITIALIZE PACKED ARRAY TO AVOID SURPRISES.

      JPAK(1:PDIM) = ZERO

!     * THIS SEQUENCE IS DONE TO FILL JPAK THE FIRST TIME AND
!     * AVOID SHIFTING BY NOTHING.

      DO  K = UN,DIM,PACK
          JJ       = (K-1)/PACK + 1
          JPAK(JJ) = J(K)
      END DO

!     * AND THEN, WE CONTINUE THE PROCESSING UP TO (UNPACK-1) TIMES
!     * BY SHIFTING JPAK TO THE LEFT AND TRANSFERING "BITS" BITS AT
!     * A TIME FROM J BY INCREMENT OF "PACK" J-ELEMENT POSITIONS.

      DO  POINT = DEUX,PACK-1
          DO  K = POINT,DIM,PACK
              JJ       = K/PACK + 1
              JPAK(JJ) = OR( ISHFT( JPAK(JJ), BITS), J(K))
          END DO
      END DO

      IF (PACK /= 1)                                           THEN
          DO  K = PACK,DIM,PACK
              JJ       = K/PACK
              JPAK(JJ) = OR( ISHFT( JPAK(JJ), BITS), J(K))
          END DO
      END IF

!     * ZEROS FILL UP THE UNOCCUPIED LOWER ORDER BITS IN JPAK LAST
!     * ELEMENT BY SHIFTING (PACK-RMNDR)*BITS COLUMNS TO THE LEFT.
!     * (I.E. LEFT JUSTIFY THE DATA BITS IN THE LAST ELEMENT OF JPAK)

      IF (RMNDR /= 0) 
     +    JPAK(PDIM) = ISHFT( JPAK(PDIM), (PACK-RMNDR)*BITS)

      RETURN

!--------------------------------------------------------------------
      END 

      SUBROUTINE sbytes1 (JPAK,X,BITS,DIM,XSCAL,XMIN)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION JPAK(*)
      REAL(4)   X(*),XSCAL,XMIN

!     * SEPTEMBER 16/1991 - B.DUGAS (RPN)

!     * THIS IS A SUBSET OF THE CRAY  SBYTESB ROUTINE (WRITTEN BY 
!     * J.F. FORTIN, CCRN - DECEMBER/1983).  IN THIS VERSION, THE 
!     * OUTPUT HAS TO FILL 32-BIT INTEGERS, IF THE RESULTS ARE TO
!     * REMAIN  IDENTICAL WITH THE CRAY VERSION.  THIS IS THE 1/1
!     * PACKING DENSITY VERSION.

!     * THIS SUBROUTINE PACKS AN ARRAY X OF  DIMENSION IDIM, CONTAINING
!     * REALS INTO AN ARRAY JPAK WHICH  DIMENSION IS <= DIM.

!----------------------------------------------------------------------
!     * "BW" IS THE NUMBER OF BITS PER INTEGER WORD.  THIS VALUE 
!     * IS MACHINE  AND IMPLEMENTATION DEPENDENT. FOR EXAMPLE, A
!     * SINGLE PRECISION NUMBER WARRANTS BW = 64 ON A CRAY COMPU-
!     * TER.

      COMMON /MACHSPEC/ BW 

! ------------------------------------------------------------------- 
!     * SHIFT JPAK TO THE LEFT AND TRANSFERING "BITS" BITS AT
!     * A TIME FROM X BY INCREMENT OF "PACK" J-ELEMENT POSITIONS.

      JPAK(1:DIM) =     NINT( XSCAL * (X(1:DIM) - XMIN) )

      RETURN

!--------------------------------------------------------------------
      END 

      SUBROUTINE sbytes2 (JPAK,X,BITS,DIM,XSCAL,XMIN)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION JPAK(*)
      REAL(4)   X(2,*),XSCAL,XMIN

!     * SEPTEMBER 16/1991 - B.DUGAS (RPN)

!     * THIS IS A SUBSET OF THE CRAY  SBYTESB ROUTINE (WRITTEN BY 
!     * J.F. FORTIN, CCRN - DECEMBER/1983).  IN THIS VERSION, THE 
!     * OUTPUT HAS TO FILL 32-BIT INTEGERS, IF THE RESULTS ARE TO
!     * REMAIN  IDENTICAL WITH THE CRAY VERSION.  THIS IS THE 2/1
!     * PACKING DENSITY VERSION.

!     * THIS SUBROUTINE PACKS AN ARRAY X OF  DIMENSION IDIM, CONTAINING
!     * REALS INTO AN ARRAY JPAK WHICH  DIMENSION IS <= DIM.

!----------------------------------------------------------------------
!     * "BW" IS THE NUMBER OF BITS PER INTEGER WORD.  THIS VALUE 
!     * IS MACHINE  AND IMPLEMENTATION DEPENDENT. FOR EXAMPLE, A
!     * SINGLE PRECISION NUMBER WARRANTS BW = 64 ON A CRAY COMPU-
!     * TER.

      COMMON /MACHSPEC/ BW 

!     * (TRANSLATING) STATEMENT FUNCTIONS.
#     if defined (HP) || defined (SUN)
      AND( NUM1, NUM2) = IAND   ( NUM1, NUM2)
      OR ( NUM1, NUM2) = IOR    ( NUM1, NUM2)
#     endif
!     ISHL( NUM,  BITS) = ISHFT (NUM,  BITS)

! ------------------------------------------------------------------- 
!     * PACKING DENSITY.

      PACK = 2
      DIM2 = DIM/PACK

!     * REMAINDER OF ELEMENTS WHICH REQUIRE BUT MAY NOT FILL A JPACK
!     * ELEMENT SPACE.

      RMNDR = MOD( DIM, PACK )

!     * SHIFT JPAK TO THE LEFT AND TRANSFERING "BITS" BITS AT
!     * A TIME FROM X BY INCREMENT OF "PACK" J-ELEMENT POSITIONS.

      DO  K = 1,DIM2

          JPAK(K) =     NINT( XSCAL * (X(1,K) - XMIN) )

          JPAK(K) = OR( 
     +                  ISHFT( JPAK(K),BITS ),
     +                  NINT( XSCAL * (X(2,K) - XMIN) ) 
     +                )

      END DO

!     * ZEROS FILL THE UNOCCUPIED LOWER ORDER BITS IN JPAK LAST ELEMENT.

      MASQUE1 = ISHFT( 2**BITS-1,1*BITS ) 

      IF (RMNDR == 1)
     +    JPAK(K) = AND( 
     +                   ISHFT( 
     +                          NINT( XSCAL * (X(1,K) - XMIN) ),
     +                          1*BITS
     +                        ), 
     +                   MASQUE1
     +                 )

      RETURN

!--------------------------------------------------------------------
      END 

      SUBROUTINE sbytes3 (JPAK,X,BITS,DIM,XSCAL,XMIN)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION JPAK(*)
      REAL(4)   X(3,*),XSCAL,XMIN

!     * SEPTEMBER 16/1991 - B.DUGAS (RPN)

!     * THIS IS A SUBSET OF THE CRAY  SBYTESB ROUTINE (WRITTEN BY 
!     * J.F. FORTIN, CCRN - DECEMBER/1983).  IN THIS VERSION, THE 
!     * OUTPUT HAS TO FILL 32-BIT INTEGERS, IF THE RESULTS ARE TO
!     * REMAIN  IDENTICAL WITH THE CRAY VERSION.  THIS IS THE 3/1
!     * PACKING DENSITY VERSION.

!     * THIS SUBROUTINE PACKS AN ARRAY X OF  DIMENSION IDIM, CONTAINING
!     * REALS INTO AN ARRAY JPAK WHICH  DIMENSION IS <= DIM.

!----------------------------------------------------------------------
!     * "BW" IS THE NUMBER OF BITS PER INTEGER WORD.  THIS VALUE 
!     * IS MACHINE  AND IMPLEMENTATION DEPENDENT. FOR EXAMPLE, A
!     * SINGLE PRECISION NUMBER WARRANTS BW = 64 ON A CRAY COMPU-
!     * TER.

      COMMON /MACHSPEC/ BW 

!     * (TRANSLATING) STATEMENT FUNCTIONS.
#     if defined (HP) || defined (SUN)
      AND( NUM1, NUM2) = IAND   ( NUM1, NUM2)
      OR ( NUM1, NUM2) = IOR    ( NUM1, NUM2)
#     endif
!     ISHL( NUM,  BITS) = ISHFT (NUM,  BITS)

! ------------------------------------------------------------------- 
!     * PACKING DENSITY.

      PACK = 3
      DIM3 = DIM/PACK

!     * REMAINDER OF ELEMENTS WHICH REQUIRE BUT MAY NOT FILL A JPACK
!     * ELEMENT SPACE.

      RMNDR = MOD( DIM, PACK )

!     * SHIFT JPAK TO THE LEFT AND TRANSFERING "BITS" BITS AT
!     * A TIME FROM X BY INCREMENT OF "PACK" J-ELEMENT POSITIONS.

      DO  K = 1,DIM3

          JPAK(K) =     NINT( XSCAL * (X(1,K) - XMIN) )

          JPAK(K) = OR( 
     +                  ISHFT( JPAK(K),BITS ),
     +                  NINT( XSCAL * (X(2,K) - XMIN) ) 
     +                )

          JPAK(K) = OR( 
     +                  ISHFT( JPAK(K),BITS ),
     +                  NINT( XSCAL * (X(3,K) - XMIN) ) 
     +                )

      END DO

!     * ZEROS FILL THE UNOCCUPIED LOWER ORDER BITS IN JPAK LAST ELEMENT.

      MASQUE1 = ISHFT( 2**(1*BITS)-1,2*BITS )
      MASQUE2 = ISHFT( 2**(2*BITS)-1,1*BITS )

      IF (RMNDR == 1)
     +    JPAK(K) = AND( 
     +                   ISHFT( 
     +                          NINT( XSCAL * (X(1,K) - XMIN) ),
     +                          2*BITS
     +                        ), 
     +                   MASQUE1
     +                 )


      IF (RMNDR == 2) 
     +    JPAK(K) = AND( 
     +                   OR( 
     +                       ISHFT(
     +                              NINT( XSCAL * (X(1,K) - XMIN) ),
     +                              2*BITS
     +                            ), 
     +                       ISHFT(
     +                              NINT( XSCAL * (X(2,K) - XMIN) ),
     +                              1*BITS
     +                            )
     +                     ),
     +                   MASQUE2
     +                 )

      RETURN

!--------------------------------------------------------------------
      END 

      SUBROUTINE sbytes4 (JPAK,X,BITS,DIM,XSCAL,XMIN)

      IMPLICIT INTEGER(4) (A-Z)

      DIMENSION JPAK(*)
      REAL(4)   X(4,*),XSCAL,XMIN

!     * SEPTEMBER 16/1991 - B.DUGAS (RPN)

!     * THIS IS A SUBSET OF THE CRAY  SBYTESB ROUTINE (WRITTEN BY 
!     * J.F. FORTIN, CCRN - DECEMBER/1983).  IN THIS VERSION, THE 
!     * OUTPUT HAS TO FILL 32-BIT INTEGERS, IF THE RESULTS ARE TO
!     * REMAIN  IDENTICAL WITH THE CRAY VERSION.  THIS IS THE 4/1
!     * PACKING DENSITY VERSION.

!     * THIS SUBROUTINE PACKS AN ARRAY X OF  DIMENSION IDIM, CONTAINING
!     * REALS INTO AN ARRAY JPAK WHICH  DIMENSION IS <= DIM.

!----------------------------------------------------------------------
!     * "BW" IS THE NUMBER OF BITS PER INTEGER WORD.  THIS VALUE 
!     * IS MACHINE  AND IMPLEMENTATION DEPENDENT. FOR EXAMPLE, A
!     * SINGLE PRECISION NUMBER WARRANTS BW = 64 ON A CRAY COMPU-
!     * TER.

      COMMON /MACHSPEC/ BW 

!     * (TRANSLATING) STATEMENT FUNCTIONS.
#     if defined (HP) || defined (SUN)
      AND( NUM1, NUM2) = IAND   ( NUM1, NUM2)
      OR ( NUM1, NUM2) = IOR    ( NUM1, NUM2)
#     endif
!     ISHL( NUM,  BITS) = ISHFT (NUM,  BITS)

! ------------------------------------------------------------------- 
!     * PACKING DENSITY.

      PACK = 4
      DIM4 = DIM/PACK

!     * REMAINDER OF ELEMENTS WHICH REQUIRE BUT MAY NOT FILL A JPACK
!     * ELEMENT SPACE.

      RMNDR = MOD( DIM, PACK )

!     * SHIFT JPAK TO THE LEFT AND TRANSFERING "BITS" BITS AT
!     * A TIME FROM X BY INCREMENT OF "PACK" J-ELEMENT POSITIONS.

      DO  K = 1,DIM4

          JPAK(K) =     NINT( XSCAL * (X(1,K) - XMIN) )

          JPAK(K) = OR( 
     +                  ISHFT( JPAK(K),BITS ),
     +                  NINT( XSCAL * (X(2,K) - XMIN) ) 
     +                )

          JPAK(K) = OR( 
     +                  ISHFT( JPAK(K),BITS ),
     +                  NINT( XSCAL * (X(3,K) - XMIN) ) 
     +                )

          JPAK(K) = OR( 
     +                  ISHFT( JPAK(K),BITS ),
     +                  NINT( XSCAL * (X(4,K) - XMIN) ) 
     +                )

      END DO

!     * ZEROS FILL THE UNOCCUPIED LOWER ORDER BITS IN JPAK LAST ELEMENT.

      MASQUE1 = ISHFT( 2**(1*BITS)-1,3*BITS )
      MASQUE2 = ISHFT( 2**(2*BITS)-1,2*BITS )
      MASQUE3 = ISHFT( 2**(3*BITS)-1,1*BITS )

      IF (RMNDR == 1)
     +    JPAK(K) = AND( 
     +                   ISHFT( 
     +                          NINT( XSCAL * (X(1,K) - XMIN) ),
     +                          3*BITS
     +                        ), 
     +                   MASQUE1
     +                 )

      IF (RMNDR == 2) 
     +    JPAK(K) = AND( 
     +                   OR( 
     +                       ISHFT(
     +                              NINT( XSCAL * (X(1,K) - XMIN) ),
     +                              3*BITS
     +                            ), 
     +                       ISHFT(
     +                              NINT( XSCAL * (X(2,K) - XMIN) ),
     +                              2*BITS
     +                            )
     +                     ),
     +                   MASQUE2
     +                 )

      IF (RMNDR == 3) 
     +    JPAK(K) = AND( 
     +                   OR(
     +                       OR( 
     +                           ISHFT( 
     +                                  NINT( XSCAL*(X(1,K) - XMIN) ),
     +                                  3*BITS
     +                                ), 
     +                           ISHFT( 
     +                                  NINT( XSCAL*(X(2,K) - XMIN) ),
     +                                  2*BITS
     +                                )
     +                         ),
     +                       ISHFT(
     +                              NINT( XSCAL*(X(3,K) - XMIN) ),
     +                              1*BITS
     +                            )
     +                     ),
     +                   MASQUE3 
     +                 )

      RETURN

!--------------------------------------------------------------------
      END 

