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
C     $Log: pacc92.ftn,v $
C     Revision 3.11  2014/12/03 23:26:35  dugas
C     Enlever les enonces EQUIVALENCE.
C
C     Revision 3.10  2014/09/25 18:42:03  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.9  2013/02/07 21:32:55  bernard
C     Ne plus generer d'erreur lorsque NBITS=1.
C
C     Revision 3.8  2008/02/06 19:43:24  dugas
C     Allocation dynamique simplifiee et autre retour en arriere.
C
C     Revision 3.6  2008/01/15 16:48:21  dugas
C     Retour a la version precedente.
C
C     Revision 3.4  1999/04/08 20:29:03  armnrbd
C     Remettre le common bloc ZZVERBO dans SBYT92.
C
C     Revision 3.3  1999/04/08 19:28:45  armnrbd
C     Utiliser le comdeck MACHTYPE.CDK.
C
C     Revision 3.2  1995/07/04 19:58:35  armnrbd
C     Allouer IXX un peu plus grand afin de s'assurer que le
C     slack (64-bits) soit bien initialise a zero dans SBYT92.
C
C     Revision 3.1  94/11/17  14:13:54  14:13:54  armnrbd (Bernard Dugas)
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C     
C     Revision 3.0  94/11/17  13:55:58  13:55:58  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:32:03  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.3  93/09/16  23:18:17  23:18:17  armnrbd (Bernard Dugas)
C     Modifier la declaration de IXX dans la routine PACC92.
C     
C     Revision 1.2  93/02/22  17:21:42  17:21:42  armnrbd (Bernard Dugas)
C     Utiliser HPALLOC/HPDEALLC plutot que MEMOIRH dans PACC92.
C     
C     Revision 1.1  92/11/15  14:57:02  armnrbd
C     Troisieme iteration PKTYP.
C     
C     Revision 1.0  92/11/11  16:03:22  armnrbd
C     Initial revision
C     
      SUBROUTINE PACC92(X,IX,IXPAK,NWDS,NBITS,XMIN,XMAX,KIND) 

***    MAY 22/92 - A.J. STACEY - RELEASED (PACCRN2 RENAMED TO PACC92).
***    MAY 15/92 - A.J. STACEY - BACK OUT INTERPOLATION BETWEEN XMIN
***                              AND XMAX.  OPTIMIZE LOOP TEST FOR XMAX
***                              CORRECTION.
***    APR 30/92 - A.J. STACEY - MODIFY PACCRN TO USE "XMIN" AND "XMAX"
***                              INSTEAD OF "XMIN" AND "XSCAL" AS THE 
***                              CONSTANTS THAT ARE STORED WITH THE 
***                              PACKED DATA.
***    APR 15/92 - A.J. STACEY - PACCRN RENAMED TO PACCRN2.
***    MAR 23/92 - A.J. STACEY - EXTENSIVE CHANGES TO ALLOW PACKER TO 
***                              RUN ON BOTH 32-BIT AND 64-BIT
***                              COMPUTERS.
***    JUN 14/89 - F. MAJAESS (CHANGE GBYTES/SBYTES CALL TO GBYT92/
***                            SBYT92) 
***    FEB 06/84 - R.LAPRISE, J.F.FORTIN.

***    FOR KIND.GE.0, PACK NWDS OF FLOATING POINT ARRAY X
***                   INTO ARRAY IXPAK.
***    FOR KIND.LT.0, UNPACK IXPAK INTO NWDS OF X. 

***    PACKING STEPS, (FOR EACH ARRAY ELEMENT) 
***       1- FIND THE MINIMUM (XMIN) AND MAXIMUM (XMAX)
***          VALUE IN THE ARRAY X. 
***       2- REMOVE BIAS (SUBTRACT XMIN FROM X). 
***       3- SCALE BY FACTOR (2**NBITS-1)/RANGE SO THAT
***          THE VALUES ARE BETWEEN 0. AND 1.*2**NBITS-1.
***       4- TRANSFORM INTO INTEGER ARRAY. 
***       5- PACK THIS ARRAY WITH SBYT92. 

***    UNPACKING, (FOR EACH ARRAY ELEMENT) 
***       1- UNPACK WITH GBYTESC.
***       2- TRANSFORM TO REAL ARRAY.
***       3- RESCALE FIELD WITH XMIN AND XSCAL.

***    PARAMETERS, 
***       X     = NON-PACKED ARRAY. (INPUT IF KIND .GE. 0, OUTPUT IF 
***                                  KIND .LT. 0). 
***       IX    = NOT USED
***       IXPAK = PACKED ARRAY. (OUTPUT IF KIND .GE. 0, INPUT IF 
***                              KIND .LT. 0)
***       NWDS  = NUMBER OF ELEMENTS IN ARRAY X (INPUT). 
***       NBITS = NUMBER OF BITS TO CODE THE ARRAY ELEMENTS IN MEMORY. 
***       XMIN  = MINIMUM OF X ARRAY. (KIND .GE. 0=OUTPUT,ELSE INPUT)
***       XSCAL = SCALE USED = (2**NBITS-1)/(XMAX-XMIN). (OUTPUT KIND
***                                               .GE. 0, INPUT ELSE)
***       KIND  = PACKING (KIND .GE. 0) OR UNPACKING (KIND .LT. 0).

      IMPLICIT      none

      INTEGER       NWDS
      REAL          X(NWDS),XMIN,XMAX
      REAL*8        XSCALI,XTMP,VLARGE,RANGE,
     +              XSCAL,SMALL,BIGGEST,BIGSTI,TEST
      INTEGER       I,IX,IXPAK(1),LARGEST,
     +              NBITS,KIND,MAXINT,IER

      INTEGER,      SAVE :: NWDS0 = -1
      INTEGER,      DIMENSION(:), ALLOCATABLE,SAVE :: IXX

#     include      "machtype.cdk" 

***    DECLARE AN ARRAY TO IMPLEMENT TRUNCATION OF FLOATING POINT NUMBS.

      INTEGER            IXSCALI     ,       IXTMP
      POINTER     ( IXSI,IXSCALI(2) ),( IXTP,IXTMP(2) )

      EXTERNAL      XIT,SBYT92,GBYT92

      DATA          MAXINT /  32   /
      DATA          SMALL,  VLARGE /
     +              1.D-38, 1.D+38 /
* ------------------------------------------------------------------- 
      ! EQUIVALENCE ( IXSCALI,   XSCALI ), ( IXTMP,   XTMP )
      IXSI = LOC( XSCALI ) ; IXTP = LOC( XTMP )

***    ARITHMETIC IS NOT GUARANTEED IF NBITS IS TOO LARGE. 

      IF (NBITS .GT. MAXINT) CALL XIT (' PACC92 ',-1 ) 
      IF (NBITS .LE. 0)      CALL XIT (' PACC92 ',-2 )

***    CALCULATE LARGEST INTEGER THAT CAN BE REPRESENTED BY "NBITS".
***    THE ORIGINAL CODE WON'T WORK ON 32-BIT MACHINE SINCE 2**NBITS CANNOT BE 
***    COMPUTED IF NBITS=32 (WHICH IS A VERY LIKELY CHOICE OF 
***    PACKING FACTOR).  NOTE ALSO THAT ON 32 BIT MACHINES, THE LARGEST
***    REPRESENTABLE *POSITIVE* INTEGER CANNOT HAVE A ONE IN BIT POSITION 31.
***    (THIS IS THE SIGN BIT!)
***    THUS FOR A PACKING FACTOR OF NPACK64=2, ONLY 31 BITS CAN BE USED
***    TO STORE PACKED DATA.  FOR ANY PACKING FACTOR GREATER THAN 2, ALL THE
***    BITS DEFINED BY NBITS MAY BE USED TO REPRESENT FLOATING POINT DATA.
***   
***    NOTE THAT FOR NPACK=2 "BIGGEST" WILL BE SLIGHTLY DIFFERENT DEPENDING WHETHER
***    FLOATING POINT NUMBERS ARE 32-BITS OR 64-BITS BECAUSE THE "FLOAT"
***    OPERATION WILL TRUNCATE SOME OF THE LEAST SIGNIFICANT BITS IN THE 
***    LARGEST 31-BIT INTEGER.  THIS WILL LEAD TO A SYSTEMATICALLY SMALLER
***    SET OF DATA IF DATA IS PACKED BY A 64-BIT MACHINE AND ANALYZED IN
***    32-BITS.  THIS WILL CAUSE AN ERROR IN (APPROX) THE 7-TH SIGNIFICANT
***    DIGIT FOR NPACK=2.  OTHER PACKING FACTORS (>2) ARE STILL OK.

      LARGEST = IBITS( -1,0,MIN( 31,NBITS ) )
      BIGGEST = DBLE( LARGEST )

      IF (NWDS0.NE.NWDS)                                       THEN

***        ALLOCATE SUFFICIENT (INTEGER) WORKING MEMORY.

          IF (NWDS0.NE.-1) 
     +    DEALLOCATE( IXX )
            ALLOCATE( IXX((NWDS+2)*INTSIZE) ) 

          NWDS0 = NWDS

      END IF

***    TEST FOR PACKING OR UNPACKING CASE. 

* ------------------------------------------------------------------- 
      IF (KIND .GE. 0)                                         THEN
* ------------------------------------------------------------------- 
***        PACKING CASE. 

***        FIND EXTREMA. 

          XMIN =  VLARGE
          XMAX = -VLARGE

          DO I=1,NWDS
              XMIN = MIN( X(I),XMIN )
              XMAX = MAX( X(I),XMAX )
          END DO

                                     RANGE = XMAX - XMIN 
          IF (ABS( RANGE ).LE.SMALL) RANGE = 1.

***        SCALE FIELD AND REPRESENT AS INTEGER. 

          XSCAL = BIGGEST / RANGE 

          IF (INTSIZE.EQ.2)                                    THEN
              DO I = 1,NWDS
                  IXX(2*I  ) = NINT( XSCAL * (X(I) - XMIN) )
              END DO
              DO I = 1,NWDS
                  IXX(2*I-1) = 0
              END DO
          ELSE IF (INTSIZE.EQ.1)                              THEN
              DO I = 1,NWDS 
                  IXX(I) = NINT( XSCAL * (X(I) - XMIN) ) 
              END DO
          END IF

***       PACK INTEGER VALUES INTO NBITS. 

          CALL SBYT92( IXPAK,IXX,NBITS,NWDS )

* ------------------------------------------------------------------- 
      ELSE
* ------------------------------------------------------------------- 
***        UNPACKING CASE. 
***        EXPAND EVERY NBITS OF IXPAK INTO WORD OF IXX.

          CALL GBYT92( IXPAK,IXX,NBITS,NWDS )

***        RESCALE X FIELD. (NOTE THAT "IXX" MAY BE A 32-BIT INTEGER ARRAY.)

                                     RANGE = XMAX - XMIN 
          IF (ABS( RANGE ).LE.SMALL) RANGE = 1.

          TEST   = RANGE/BIGGEST/2.0
          XSCALI = RANGE / BIGGEST
          BIGSTI = 1./BIGGEST

          IF (INTSIZE .EQ. 1)                                   THEN

              DO 300 I = 1,NWDS 

                  X(I) = DBLE( IXX(I) ) * XSCALI + XMIN

***                THE FOLLOWING CODE TWEAKS THE VALUE OF X(I) TO EXACTLY 
***                EQUAL XMAX WHEN X(I)=XMAX WITHIN ROUNDOFF ERROR.  THIS
***                CORRECTS A SINGLE BIT ERROR ON THE SX-3 DUE TO ROUNDOFF 
***                ERROR. IT IS NOT NECESSARY TO PERFORM THE SAME TWEAKING 
***                FOR "XMIN" IF "X(I) = FLOAT(IX(I))*XSCALI + XMIN" IS USED.

                  IF (ABS( X(I)-XMAX ).LT.TEST) X(I) = XMAX

  300         CONTINUE

          ELSE IF (INTSIZE.EQ.2)                                   THEN

              DO 310 I = 1,NWDS 

                  X(I) = DBLE( IXX(2*I) ) * XSCALI + XMIN
                  IF (ABS( X(I)-XMAX ).LT.TEST) X(I) = XMAX

  310         CONTINUE

          END IF

* ------------------------------------------------------------------- 
      END IF

      RETURN

* ------------------------------------------------------------------- 
      END 

      SUBROUTINE GBYT92(JPAK,J,NBITS,IDIM)

***    JUN 16/92 - A.J. STACEY - MODIFIED TO HANDLE PARTIALLY PACKED WORDS
***                              WITHOUT OVERWRITING.
***    APR 15/92 - A.J. STACEY - REWRITTEN GBYTESC ROUTINE UNPACKS
***                              INTEGERS.
***   
***    PURPOSE: GBYT92 UNPACKS INTEGERS FROM EITHER 32- OR 64-BIT
***             WORDS.
***             THIS ROUTINE IS COMPLETELY PORTABLE BETWEEN 32
***             AND 64-BIT MACHINES.
***   
***             THIS ROUTINE IS THE INVERSE OF SBYT92.
***    

      IMPLICIT  INTEGER (A-Z)

      INTEGER   J(1),JPAK(1)

      LOGICAL            INFO
      COMMON   /ZZVERBO/ INFO

#     include  "machtype.cdk" 

*-------------------------------------------------------------------------------
      NPACK32   = 32/NBITS
      MASK      = IBITS(-1,0,MIN(31,NBITS))
      IF (NPACK32 .LE. 0) THEN
         IF (INFO) WRITE(6,6100)
 6100    FORMAT('0*** ERROR *** GBYT92: NBITS > 32 IS ILLEGAL.')
         CALL XIT(' Gbyt92 ',-1 )
      END IF
      NPACK64        = 2*NPACK32

***    THIS IS THE UNOPTIMIZED CODE.

*     JJ               = 1
*     II               = 1
*     DO 500 K=1,LEN64*MACHINE
*       DO 400 M=2/MACHINE-1,0,-1
*         DO 400 I=NPACK32-1,0,-1
*           IOFFSET       = M*32
*           J(II*INTSIZE) = IBITS(JPAK(JJ),I*NBITS+IOFFSET,NBITS)
*           II            = II + 1
* 400   CONTINUE
*       JJ             = JJ + 1
* 500 CONTINUE

***    THIS IS THE OPTIMIZED CODE.  NOTE THE LOOP INVERSION.

      IRMDR1      = MOD(IDIM,NPACK64)
      IRMDR2      = MOD(IRMDR1,NPACK32)
      LEN64       = IDIM/NPACK64
      MLAST       = IRMDR1/NPACK32
      II          = 1
      INC         = NPACK32*2/MACHINE
      IF (MACHINE .EQ. 1) THEN            ! 64-BIT MACHINE.

        DO 400 M=1,0,-1
        DO 400 I=NPACK32-1,0,-1
          IPT = I*NBITS + M*32
          KPT = II
          DO 500 K=1,LEN64
*           J(KPT) = IBITS(JPAK(K),IPT,NBITS)!FAILS ON MIPS (SIGN BIT RETAINED)
            J(KPT) = IAND(ISHFT(JPAK(K),-IPT),MASK)
            KPT    = KPT + INC
  500     CONTINUE
          II = II + 1
  400   CONTINUE

        K   = LEN64 + 1
        IF (IRMDR1 .NE. 0) THEN            ! UNPACK PARTIAL WORD OF DATA.
   	  KPT = LEN64*NPACK64 + 1
          IF (MLAST .GT. 0) THEN           ! UNPACK NPACK32 PARCELS.
            DO 520 I=NPACK32-1,0,-1
             IPT    = I*NBITS + 32
*            J(KPT) = IBITS(JPAK(K),IPT,NBITS)!FAILS ON MIPS (SIGN BIT RETAINED)
             J(KPT) = IAND(ISHFT(JPAK(K),-IPT),MASK)
             KPT    = KPT + 1
  520       CONTINUE
          END IF

***        UNPACK REMAINDER OF  DATA.

          MLEFT = 0
          IF (MLAST .EQ. 0) MLEFT = 1
          DO 540 I=NPACK32-1,NPACK32-IRMDR2,-1
            IPT    = I*NBITS + MLEFT*32
*           J(KPT) = IBITS(JPAK(K),IPT,NBITS)!FAILS ON MIPS (SIGN BIT RETAINED)
            J(KPT) = IAND(ISHFT(JPAK(K),-IPT),MASK)
            KPT    = KPT + 1
  540     CONTINUE
        END IF

      ELSE IF (MACHINE .EQ. 2) THEN        ! 32-BIT MACHINE.

        LEN32  = IDIM/NPACK32
        IRMDR2 = MOD(IDIM,NPACK32)
        INC2   = INC*INTSIZE
        DO 700 I=NPACK32-1,0,-1
          IPT  = I*NBITS
          KPT  = II*INTSIZE
          DO 800 K=1,LEN32
*           J(KPT) = IBITS(JPAK(K),IPT,NBITS)!FAILS ON MIPS (SIGN BIT RETAINED)
            J(KPT) = IAND(ISHFT(JPAK(K),-IPT),MASK)
            KPT    = KPT + INC2
  800     CONTINUE
          II = II + 1
  700   CONTINUE

        IF (IRMDR2 .NE. 0) THEN          ! UNPACK REMAINDER OF  DATA.
          K        = LEN32 + 1
          KPT      = (IDIM - IRMDR2 + 1)*INTSIZE
          DO 840 I=NPACK32-1,NPACK32-IRMDR2,-1
            IPT    = I*NBITS
*           J(KPT) = IBITS(JPAK(K),IPT,NBITS) !FAILS ON MIPS (SIGN BIT RETAINED)
            J(KPT) = IAND(ISHFT(JPAK(K),-IPT),MASK)
            KPT    = KPT + INTSIZE
  840     CONTINUE
        END IF

      ELSE

         IF (INFO) WRITE(6,6200)
 6200    FORMAT('0*** ERROR *** GBYT92: ILLEGAL MACHINE TYPE.')
         CALL XIT(' Gbyt92 ',-2 )

      END IF

      RETURN
      END

      SUBROUTINE SBYT92(JPAK,J,NBITS,IDIM)

***    APR 15/92 - A.J. STACEY - REWRITTEN SBYTESB ROUTINE PACKS
***                              INTEGERS.

***    PURPOSE: SBYT92 PACKS INTEGERS INTO EITHER 32- OR 64-BIT WORDS
***             SUCH THAT NO PACKED DATA ELEMENT SPANS EITHER A 32 OR 
***             64 BIT BOUNDARY. THIS ROUTINE IS COMPLETELY PORTABLE
***             BETWEEN 32 AND 64-BIT MACHINES.

***    NOTES:
***    1) IDIM IS THE NUMBER OF FLOATING POINT NUMBERS THAT
***       ARE BEING PACKED.  THE INTEGER EQUIVALENTS ARE BEING 
***       STORED IN J, WHICH CAN BE EITHER A 32-BIT OR 64-BIT INTEGER
***       ARRAY DEPENDING ON WHETHER THE NATIVE MACHINE IS A 32-BIT 
***       MACHINE OR A 64-BIT MACHINE.
***    2) NO CHECK ON INTEGERS IN J BEING VALID INTEGERS IN RANGE 
***       (0,2**NBITS-1)
***    3) IF THE NUMBER OF WORDS IN THE "J" ARRAY IS NOT DIVISIBLE
***       BY THE PACKING FACTOR, THEN THE LAST WORD OF THE "JPAK"
***       ARRAY CONTAINS PARTIAL DATA.
***    4) WE ADOPT THE CONVENTION THAT THE INPUT ARRAY "J" CONTAINS
***       SIGNIFICANT DATA ONLY IN THE LOW-ORDER 32 BITS OF EACH 64-BIT
***       FIELD.  ON 32-BIT MACHINES THIS MEANS THAT EACH EVEN ELEMENT
***       OF "J" CONTAINS THE DATA TO BE PACKED.

      IMPLICIT  INTEGER (A-Z)

      INTEGER   J(1),JPAK(1)

      LOGICAL            INFO
      COMMON   /ZZVERBO/ INFO 

#     include  "machtype.cdk" 

*-------------------------------------------------------------------------------
      NPACK32   = 32/NBITS
      MASK      = IBITS(-1,0,MIN(31,NBITS))
      IF (NPACK32 .LE. 0) THEN
         IF (INFO) WRITE(6,6100)
 6100    FORMAT('0*** ERROR *** SBYT92: NBITS > 32 IS ILLEGAL.')
         CALL XIT(' Sbyt92 ',-1 )
      END IF
      NPACK64          = 2*NPACK32
      IRMNDR           = MOD(IDIM,NPACK64)
      IF (IRMNDR .EQ. 0) THEN
          JPDIM        = IDIM/NPACK64
      ELSE
          JPDIM        = (IDIM/NPACK64) + 1
      END IF
      DO 100 K=1,JPDIM*MACHINE
        JPAK(K)        = 0
  100 CONTINUE
      LEN64          = ((IDIM/NPACK64) +
     1  	       ((MOD(IDIM,NPACK64)+NPACK64-1)/NPACK64))

***    THIS IS THE UNOPTIMIZED CODE.

*     JJ             = 1
*     II             = 1
*     DO 500 K=1,LEN64*MACHINE
*       DO 400 M=2/MACHINE-1,0,-1
*         DO 400 I=NPACK32-1,0,-1
*           IOFFSET  = M*32
*           JPAK(JJ) = IOR(
*    1  		   JPAK(JJ),
*    2  		   ISHFT(J(II*INTSIZE),I*NBITS+IOFFSET)
*    3  		  )
*           II       = II + 1
* 400   CONTINUE
*       JJ           = JJ + 1
* 500 CONTINUE

***    THIS IS THE OPTIMIZED CODE.

      II = 1
      INC = NPACK32*2/MACHINE
      IF (MACHINE .EQ. 1) THEN         ! 64-BIT MACHINE.
      DO 400 M=1,0,-1
        DO 400 I=NPACK32-1,0,-1
          IPT = I*NBITS + M*32
          KPT = II
          DO 500 K=1,LEN64
            JPAK(K) = IOR(
     1  		   JPAK(K),
     2  		   ISHFT(IAND(J(KPT),MASK),IPT)
     3  		 )
            KPT     = KPT + INC
  500     CONTINUE
          II = II + 1
  400 CONTINUE
      ELSE IF (MACHINE .EQ. 2 ) THEN   ! 32-BIT MACHINE.
      LEN32 = LEN64*2
      INC2  = INC*INTSIZE
      DO 700 I=NPACK32-1,0,-1
          IPT = I*NBITS
          KPT = II*INTSIZE
          DO 800 K=1,LEN32
            JPAK(K) = IOR(
     1  		   JPAK(K),
     2  		   ISHFT(IAND(J(KPT),MASK),IPT)
     3  		 )
            KPT     = KPT + INC2
  800     CONTINUE
          II = II + 1
  700 CONTINUE
      ELSE
         IF (INFO) WRITE(6,6200)
 6200    FORMAT('0*** ERROR *** SBYT92: ILLEGAL MACHINE TYPE.')
         CALL XIT(' Sbyt92 ',-2 )
      END IF

      RETURN
      END

