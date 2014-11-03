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
C     $Log: pfa.ftn,v $
C     Revision 3.4  2014/09/25 18:42:03  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.3  2008/04/28 21:38:53  dugas
C     Corriger l'usage des macros pour r.gppf (passe 2).
C
C     Revision 3.2  2003/09/15 16:19:46  dugas
C     Re-definir le macro lot_maximum de 1024 a 4 suite a
C      une modification equivalente dans les versions libpriv.a
C      de RPASSM8 et QPASSM8 sur nos frontaux survenues recamment.
C
C     Revision 3.1  1994/11/17 14:13:59  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:56:02  13:56:02  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:32:06  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.1  93/02/16  10:58:56  armnrbd
C     Renommer setpfa,pfa,qpfa et rpfa en ...
C      ...     setpfa2,pfa2,qpfa2,rpfa2.
C     
C     Revision 1.0  92/02/21  11:34:16  armnrbd
C     Initial revision
C     
#     if !defined (lot_maximum)
#         define   lot_maximum 4
#     endif
************************************************************************
*                                                                      *
*     SUBROUTINE SETPFA2 - SET UP TABLES READY FOR PRIME FACTOR FFT    *
*                                                                      *
*     ENTRY: CALL SETPFA2(IFAX,ILIST,N,INCR,INCS,IFORM,IERR)           *
*                                                                      *
*     ON EXIT:                                                         *
*     --------                                                         *
*     IFAX(1) CONTAINS NFAX, THE NUMBER OF FACTORS OF N                *
*     IFAX(2) TO IFAX(NFAX+1) CONTAIN THE FACTORS OF N                 *
*                                                                      *
*     ILIST(1) CONTAINS (2**24)*NFOLD(1) + NFOLD(2)                    *
*     ILIST(2) CONTAINS (2**24)*NFOLD(3) + NFOLD(4)                    *
*       WHERE NFOLD(M), M=1,2,3,4, IS THE NUMBER OF M-DIMENSIONAL      *
*       "FOLDS" IN THE REAL/HALF-COMPLEX VERSION OF THE ALGORITHM      *
*     ILIST(3) TO ILIST(2*N) CONTAIN A LIST OF ADDRESSES USED          *
*       DURING THE "FOLDING" PHASE OF THE REAL/HALF-COMPLEX            *
*       VERSION OF THE ALGORITHM                                       *
*     IERR  = 0, IF ALL WENT WELL                                      *
*             1, IF ILLEGAL FACTORS WERE FOUNDY                        *
*             2, IF N CONTAINS TOO MANY FACTORS                        *
*             3, IF THE FACTORS ARE NOT MUTUALLY PRIME                 *
*                                                                      *
*                                                                      *
*     ON ENTRY:                                                        *
*     ---------                                                        *
*     N = LENGTH OF TRANSFORMS                                         *
*     INCR = INCREMENT BETWEEN SUCCESSIVE REAL-SPACE DATA VALUES       *
*            WITHIN EACH TRANSFORM                                     *
*     INCS = INCREMENT BETWEEN SUCCESSIVE FOURIER-SPACE COEFF-         *
*            ICIENTS WITHIN EACH TRANSFORM                             *
*     IFORM = +1 FOR CONVENTIONAL ORDERING OF FOURIER COEFFICIENTS     *
*           = -1 FOR "FAST-IN-PLACE" ORDERING                          *
*           =  0 FOR A COMPLEX/COMPLEX TRANSFORM, IN WHICH CASE        *
*                THE ARRAY ILIST IS NOT REQUIRED                       *
*                                                                      *
*     AUTHOR:                                                          *
*     -------                                                          *
*     CLIVE TEMPERTON    R P N    JULY 1985                            *
*                                                                      *
*     MODIFIED BY BERNARD DUGAS DECEMBER 1992 -- ADDED IERR ERROR      *
*     COMPLETION CODE OUTPUT.                                          *
*                                                                      *
************************************************************************
*
      SUBROUTINE SETPFA2(IFAX,ILIST,N,INCR,INCS,IFORM,IERR )
*
      IMPLICIT INTEGER (A-Z)
*
      INTEGER  N,INCR,INCS,IFORM,IERR
      INTEGER  IFAX(5),ILIST(2*N)
      INTEGER  ITEST(8),MM(4),NFOLD(4) 
      DATA     ITEST/16,8,4,2,9,3,5,7/
      IERR=0
*
*     FIND FACTORS OF N
*     -----------------
      NFAX=0
      J=2 
      NN=N
      DO 10 I=1,8
      IFAC=ITEST(I) 
      IF (MOD(NN,IFAC).EQ.0) THEN
        IFAX(J)=IFAC
        J=J+1
        NFAX=NFAX+1 
        NN=NN/IFAC
      ENDIF
   10 CONTINUE
*
      IFAX(1)=NFAX
      IF (NN.NE.1) THEN
        WRITE(6,20) N
        IERR=1
        RETURN
      ENDIF
   20 FORMAT('0N =',I5,5X,'CONTAINS ILLEGAL FACTOR.')
      IF (NFAX.GT.4) THEN
        WRITE(6,30) N
        IERR=2
        RETURN
      ENDIF
   30 FORMAT('0N =',I5,5X,'CONTAINS TOO MANY FACTORS.')
*
*     GENERATE LIST OF VALUES OF MM AND 
*     CHECK THAT FACTORS ARE MUTUALLY PRIME
*     ------------------------------------
*
      DO 140 K = 1 , NFAX
      IFAC = IFAX(K+1)
      M = N / IFAC
      DO 110 J = 1 , IFAC
      MJ = J * M
      IF (MOD(MJ,IFAC).EQ.1) GO TO 130
  110 CONTINUE
*
      WRITE(6,120) N,(IFAX(I+1),I=1,NFAX)
  120 FORMAT('0N =',I5,5X,'IFAX =',4I4,5X,
     +    'FACTORS ARE NOT MUTUALLY PRIME.')
      IERR=3
      RETURN
*
  130 CONTINUE
      MM(K) = MJ
  140 CONTINUE

      WRITE(6,145) IFAX(1),(IFAX(I),I=2,IFAX(1)+1)
  145 FORMAT('      PFA  NFAX=',I4,', factors ...',5I4)
*
      IF (IFORM.EQ.0) RETURN
*
*     GENERATE ILIST
*     --------------
      DO 150 I = 1 , 4
      NFOLD(I) = 0
  150 CONTINUE
      KK = 3
      NH = N / 2
*
*     LIST FOR LINES
*     --------------
      IQ = 0
      DO 230 K = 1 , NFAX
      IFAC = IFAX(K+1)
      IF (IFAC.EQ.2) GO TO 230
      IZZ = ( IFAC - 1 ) / 2
      LL = 1 - MOD(N,2) + MOD(IFAC,2)
      IA = 0
      DO 220 L = 1 , LL
      DO 210 I = 1 , IZZ
      IA = IA + MM(K)
      IF (IA.GT.N) IA = IA - N
      IB = N - IA
      JA = MIN0 (IA,IB)
      JB = N - JA
      IF (IFORM.EQ.+1) THEN
        JA = 2 * JA 
        JB = JA + 1 
      ENDIF
      ILIST(KK) = INCR * IA
      ILIST(KK+1) = INCS * JA 
      ILIST(KK+2) = INCR * IB 
      ILIST(KK+3) = INCS * JB 
      KK = KK + 4
  210 CONTINUE
      IQ = IQ + IZZ 
      IA = NH
  220 CONTINUE
  230 CONTINUE
*
      NFOLD(1) = IQ 
      IF (NFAX.EQ.1) GO TO 600
*
*     LIST FOR PLANES
*     ---------------
      IQ = 0
      DO 350 K1 = 1 , NFAX - 1
      IFAC1 = IFAX(K1+1)
      IF (IFAC1.EQ.2) GO TO 350
      IZZ = ( IFAC1 - 1 ) / 2 
      DO 340 K2 = K1 + 1 , NFAX
      IFAC2 = IFAX(K2+1)
      IF (IFAC2.EQ.2) GO TO 340
      JZZ = ( IFAC2 - 1 ) / 2 
      LL = 1 - MOD(N,2) + MOD(IFAC1*IFAC2 , 2)
      IX = 0
      DO 330 L = 1 , LL
      DO 320 J = 1 , JZZ
      IX = IX + MM(K2)
      IF (IX.GT.N) IX = IX - N
      IA = IX
      IB = IX
      DO 310 I = 1 , IZZ
      IA = IA + MM(K1)
      IF (IA.GT.N) IA = IA - N
      IB = IB - MM(K1)
      IF (IB.LT.0) IB = IB + N
      IC = N - IB
      ID = N - IA
      JA = MIN0 (IA,ID)
      JB = MIN0 (IB,IC)
      JC = N - JB
      JD = N - JA
      IF (IFORM.EQ.+1) THEN
        JA = 2 * JA 
        JB = 2 * JB 
        JC = JB + 1 
        JD = JA + 1 
      ENDIF
      ILIST(KK) = INCR * IA
      ILIST(KK+1) = INCS * JA 
      ILIST(KK+2) = INCR * IB 
      ILIST(KK+3) = INCS * JB 
      ILIST(KK+4) = INCR * IC 
      ILIST(KK+5) = INCS * JC 
      ILIST(KK+6) = INCR * ID 
      ILIST(KK+7) = INCS * JD 
      KK = KK + 8
  310 CONTINUE
      IQ = IQ + IZZ 
  320 CONTINUE
      IX = NH
  330 CONTINUE
  340 CONTINUE
  350 CONTINUE
*
      NFOLD(2) = IQ 
      IF (NFAX.EQ.2) GO TO 600
*
*     LIST FOR CUBOIDS
*     ----------------
      IQ = 0
      DO 470 K1 = 1 , NFAX - 2
      IFAC1 = IFAX(K1+1)
      IF (IFAC1.EQ.2) GO TO 470
      IZZ = ( IFAC1 - 1 ) / 2 
      DO 460 K2 = K1 + 1 , NFAX - 1
      IFAC2 = IFAX(K2+1)
      IF (IFAC2.EQ.2) GO TO 460
      JZZ = (IFAC2 - 1) / 2
      DO 450 K3 = K2 + 1 , NFAX
      IFAC3 = IFAX(K3+1)
      IF (IFAC3.EQ.2) GO TO 450
      KZZ = (IFAC3 - 1) / 2
      LL = 1 - MOD(N,2) + MOD(IFAC1*IFAC2*IFAC3 , 2)
      IX = 0
      DO 440 L = 1 , LL
      DO 430 K = 1 , KZZ
      IX = IX + MM(K3)
      IF (IX.GT.N) IX = IX - N
      IY = IX
      IZ = IX
      DO 420 J = 1 , JZZ
      IY = IY + MM(K2)
      IF (IY.GT.N) IY = IY - N
      IZ = IZ - MM(K2)
      IF (IZ.LT.0) IZ = IZ + N
      IA = IY
      IB = IY
      IC = IZ
      ID = IZ
      DO 410 I = 1 , IZZ
      IA = IA + MM(K1)
      IF (IA.GT.N) IA = IA - N
      IB = IB - MM(K1)
      IF (IB.LT.0) IB = IB + N
      IC = IC + MM(K1)
      IF (IC.GT.N) IC = IC - N
      ID = ID - MM(K1)
      IF (ID.LT.0) ID = ID + N
      IE = N - ID
      IF = N - IC
      IG = N - IB
      IH = N - IA
      JA = MIN0 (IA,IH)
      JB = MIN0 (IB,IG)
      JC = MIN0 (IC,IF)
      JD = MIN0 (ID,IE)
      JE = N - JD
      JF = N - JC
      JG = N - JB
      JH = N - JA
      IF (IFORM.EQ.+1) THEN
        JA = 2 * JA 
        JB = 2 * JB 
        JC = 2 * JC 
        JD = 2 * JD 
        JE = JD + 1 
        JF = JC + 1 
        JG = JB + 1 
        JH = JA + 1 
      ENDIF
      ILIST(KK) = INCR * IA
      ILIST(KK+1) = INCS * JA 
      ILIST(KK+2) = INCR * IB 
      ILIST(KK+3) = INCS * JB 
      ILIST(KK+4) = INCR * IC 
      ILIST(KK+5) = INCS * JC 
      ILIST(KK+6) = INCR * ID 
      ILIST(KK+7) = INCS * JD 
      ILIST(KK+8) = INCR * IE 
      ILIST(KK+9) = INCS * JE 
      ILIST(KK+10) = INCR * IF
      ILIST(KK+11) = INCS * JF
      ILIST(KK+12) = INCR * IG
      ILIST(KK+13) = INCS * JG
      ILIST(KK+14) = INCR * IH
      ILIST(KK+15) = INCS * JH
      KK = KK + 16
  410 CONTINUE
      IQ = IQ + IZZ 
  420 CONTINUE
  430 CONTINUE
      IX = NH
  440 CONTINUE
  450 CONTINUE
  460 CONTINUE
  470 CONTINUE
*
      NFOLD(3) = IQ 
      IF (NFAX.EQ.3) GO TO 600
*
*     LIST FOR HYPERCUBOIDS
*     ---------------------
      IQ = 0
      IFAC1 = IFAX(2)
      IF (IFAC1.EQ.2) GO TO 600
      IZZ = ( IFAC1 - 1 ) / 2 
      IFAC2 = IFAX(3)
      IF (IFAC2.EQ.2) GO TO 600
      JZZ = ( IFAC2 - 1 ) / 2 
      IFAC3 = IFAX(4)
      IF (IFAC3.EQ.2) GO TO 600
      KZZ = ( IFAC3 - 1 ) / 2 
      IFAC4 = IFAX(5)
      IF (IFAC4.EQ.2) GO TO 600
      LZZ = ( IFAC4 - 1 ) / 2 
      IT = 0
      DO 540 L = 1 , LZZ
      IT = IT + MM(4)
      IF (IT.GT.N) IT = IT - N
      IU = IT
      IV = IT
      DO 530 K = 1 , KZZ
      IU = IU + MM(3)
      IF (IU.GT.N) IU = IU - N
      IV = IV - MM(3)
      IF (IV.LT.0) IV = IV + N
      IW = IU
      IX = IU
      IY = IV
      IZ = IV
      DO 520 J = 1 , JZZ
      IW = IW + MM(2)
      IF (IW.GT.N) IW = IW - N
      IX = IX - MM(2)
      IF (IX.LT.0) IX = IX + N
      IY = IY + MM(2)
      IF (IY.GT.N) IY = IY - N
      IZ = IZ - MM(2)
      IF (IZ.LT.0) IZ = IZ + N
      IA = IW
      IB = IW
      IC = IX
      ID = IX
      IE = IY
      IF = IY
      IG = IZ
      IH = IZ
      DO 510 I = 1 , IZZ
      IA = IA + MM(1)
      IF (IA.GT.N) IA = IA - N
      IB = IB - MM(1)
      IF (IB.LT.0) IB = IB + N
      IC = IC + MM(1)
      IF (IC.GT.N) IC = IC - N
      ID = ID - MM(1)
      IF (ID.LT.0) ID = ID + N
      IE = IE + MM(1)
      IF (IE.GT.N) IE = IE - N
      IF = IF - MM(1)
      IF (IF.LT.0) IF = IF + N
      IG = IG + MM(1)
      IF (IG.GT.N) IG = IG - N
      IH = IH - MM(1)
      IF (IH.LT.0) IH = IH + N
      II = N - IH
      IJ = N - IG
      IK = N - IF
      IL = N - IE
      IM = N - ID
      IN = N - IC
      IO = N - IB
      IP = N - IA
      JA = MIN0 (IA,IP)
      JB = MIN0 (IB,IO)
      JC = MIN0 (IC,IN)
      JD = MIN0 (ID,IM)
      JE = MIN0 (IE,IL)
      JF = MIN0 (IF,IK)
      JG = MIN0 (IG,IJ)
      JH = MIN0 (IH,II)
      JI = N - JH
      JJ = N - JG
      JK = N - JF
      JL = N - JE
      JM = N - JD
      JN = N - JC
      JO = N - JB
      JP = N - JA
      IF (IFORM.EQ.+1) THEN
        JA = 2 * JA 
        JB = 2 * JB 
        JC = 2 * JC 
        JD = 2 * JD 
        JE = 2 * JE 
        JF = 2 * JF 
        JG = 2 * JG 
        JH = 2 * JH 
        JI = JH + 1 
        JJ = JG + 1 
        JK = JF + 1 
        JL = JE + 1 
        JM = JD + 1 
        JN = JC + 1 
        JO = JB + 1 
        JP = JA + 1 
      ENDIF
      ILIST(KK) = INCR * IA
      ILIST(KK+1) = INCS * JA 
      ILIST(KK+2) = INCR * IB 
      ILIST(KK+3) = INCS * JB 
      ILIST(KK+4) = INCR * IC 
      ILIST(KK+5) = INCS * JC 
      ILIST(KK+6) = INCR * ID 
      ILIST(KK+7) = INCS * JD 
      ILIST(KK+8) = INCR * IE 
      ILIST(KK+9) = INCS * JE 
      ILIST(KK+10) = INCR * IF
      ILIST(KK+11) = INCS * JF
      ILIST(KK+12) = INCR * IG
      ILIST(KK+13) = INCS * JG
      ILIST(KK+14) = INCR * IH
      ILIST(KK+15) = INCS * JH
      ILIST(KK+16) = INCR * II
      ILIST(KK+17) = INCS * JI
      ILIST(KK+18) = INCR * IJ
      ILIST(KK+19) = INCS * JJ
      ILIST(KK+20) = INCR * IK
      ILIST(KK+21) = INCS * JK
      ILIST(KK+22) = INCR * IL
      ILIST(KK+23) = INCS * JL
      ILIST(KK+24) = INCR * IM
      ILIST(KK+25) = INCS * JM
      ILIST(KK+26) = INCR * IN
      ILIST(KK+27) = INCS * JN
      ILIST(KK+28) = INCR * IO
      ILIST(KK+29) = INCS * JO
      ILIST(KK+30) = INCR * IP
      ILIST(KK+31) = INCS * JP
      KK = KK + 32
  510 CONTINUE
      IQ = IQ + IZZ 
  520 CONTINUE
  530 CONTINUE
  540 CONTINUE
*
      NFOLD(4) = IQ 
*
  600 CONTINUE
*
      ISHIFT = 2**24
      ILIST(1) = ISHIFT * NFOLD(1) + NFOLD(2)
      ILIST(2) = ISHIFT * NFOLD(3) + NFOLD(4)
*
      RETURN
      END 
*DECK FORCPFA
*        SUBROUTINE 'PFA2'
*        SELF-SORTING IN-PLACE PRIME FACTOR (COMPLEX) FFT
*
*        CALL PFA2(A,B,WORK,IFAX,INC,JUMP,N,LOT,ISIGN,IERR)
*
*        A IS FIRST REAL INPUT/OUTPUT VECTOR
*        B IS FIRST IMAGINARY INPUT/OUTPUT VECTOR 
*        WORK IS NOT USED
*        IFAX IS A LIST OF FACTORS OF N 
*              IFAX(1) CONTAINS THE NUMBER OF FACTORS (NFAX)
*              THE FACTORS FOLLOW IN IFAX(2) TO IFAX(NFAX+1)
*        INC IS THE INCREMENT WITHIN EACH DATA VECTOR
*        JUMP IS THE INCREMENT BETWEEN DATA VECTORS
*        N IS THE LENGTH OF THE TRANSFORMS
*        LOT IS THE NUMBER OF TRANSFORMS
*        ISIGN = +1 FOR FORWARD TRANSFORM
*                -1 FOR INVERSE TRANSFORM
*        IERR IS AN ERROR INDICATOR: 
*                  0 - TRANSFORM COMPLETED WITHOUT ERROR
*                  1 - FACTORS ARE NOT MUTUALLY PRIME
*                  2 - ILLEGAL FACTOR
*
*        FORTRAN VERSION BY CLIVE TEMPERTON RPN NOVEMBER 1986
*
*        MODIFIED BY CLIVE TEMPERTON MAY 1987 -- REDUNDANT PARAMETER
*        "WORK" ADDED TO CALLING LIST FOR COMPATIBILITY WITH CAL VERSION
*
*        MODIFIED BY BERNARD DUGAS DECEMBER 1992 -- 1) USE  AN 
*        "IMPLICIT NONE" FORMULATION, 2) USE VECTOR TEMPORARIES WITHIN
*        LOOPS, AND 3) USE SAME NORMALIZATIONS AS OTHER TRANSFORMS 
*        ROUTINES (SUCH AS FFT772)
*
*----------------------------------------------------------------------
*
*        DEFINITION OF TRANSFORM
*        -----------------------
*
*        IF (ISIGN.EQ.+1)
*        X(J) =       SUM(K=0,...,N-1)(C(K)*EXP(ISIGN*2*I*J*K*PI/N))
*        IF (ISIGN.EQ.-1)
*        X(J) = (1/N) SUM(K=0,...,N-1)(C(K)*EXP(ISIGN*2*I*J*K*PI/N))
*
*---------------------------------------------------------------------
*
*        FOR A MATHEMATICAL DEVELOPMENT OF THE ALGORITHM USED,
*        SEE: 
*
*        C TEMPERTON : "IMPLEMENTATION OF A SELF-SORTING IN-PLACE
*          PRIME FACTOR FFT ALGORITHM", 
*        JOURNAL OF COMPUTATIONAL PHYSICS VOL. 58 NO. 3 (1985), 283-299
*
*-------------------------------------------------------------------- 
C
      SUBROUTINE PFA2(A,B,WORK,IFAX,INC,JUMP,N,LOT,ISIGN,IERR)
*
      IMPLICIT NONE
*
      INTEGER  INC,JUMP,N,LOT,ISIGN,IERR
      REAL     A(*),B(*),WORK
      INTEGER  IFAX(*)
*
      REAL     AA,AB,AC,AD,AE,AF,AG,AH,AI,AJ,AK,AL,AM,AN,AO,AP
      REAL     BA,BB,BC,BD,BE,BF,BG,BH,BI,BJ,BK,BL,BM,BN,BO,BP
      REAL     T1,T2,T3,T4,T5,T6,TWOPI,ANGLE,FACT
      REAL     Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9
      INTEGER  NFAX,IFAC,I,J,K,L,M,MU,MM,MMU,NN,IX
      INTEGER  IA,IB,IC,ID,IE,IF,IG,IH,II,IJ,IK,IL,IM,IN,IO,IP
*
      REAL     SIN60,QRT5,SIN36,SIN45,SIN72
      DATA     SIN60/0.866025403784437/
      DATA     QRT5 /0.559016994374947/
      DATA     SIN36/0.587785252292473/,SIN72/0.951056516295154/
      DATA     SIN45/0.707106781186548/
C
      NFAX = IFAX(1)
*
      DO 1000 K=1,NFAX
C
      IFAC=IFAX(K+1)
      M=N/IFAC
      DO 100 J=1,IFAC
      MU=J
      MM=J*M
      IF (MOD(MM,IFAC).EQ.1) GO TO 110
  100 CONTINUE
*
*     ERROR EXIT - FACTORS NOT MUTUALLY PRIME
*     ----------
      IERR=1
      GO TO 2000
*
  110 CONTINUE
*
      FACT=1.0
      IF (ISIGN.EQ.-1) THEN
        MU = IFAC - MU
        IF (K.EQ.1) FACT=1.0/FLOAT(N)
      ENDIF
      MM = MM * INC 
      NN = N * INC
C
C     NOW COMPUTE THE ADDRESSES IA,IB ETC. AND SELECT THE
C     CODING FOR THE CURRENT FACTOR
C
      IA=1
      IB=IA+MM
      IF (IFAC.EQ.2) GO TO 200
      IC=IB+MM
      IF (IC.GT.NN) IC=IC-NN
      IF (IFAC.EQ.3) GO TO 300
      ID=IC+MM
      IF (ID.GT.NN) ID=ID-NN
      IF (IFAC.EQ.4) GO TO 400
      IE=ID+MM
      IF (IE.GT.NN) IE=IE-NN
      IF (IFAC.EQ.5) GO TO 500
      IF=IE+MM
      IF (IF.GT.NN) IF=IF-NN
      IG=IF+MM
      IF (IG.GT.NN) IG=IG-NN
      IF (IFAC.EQ.7) GO TO 700
      IH=IG+MM
      IF (IH.GT.NN) IH=IH-NN
      IF (IFAC.EQ.8) GO TO 800
      II=IH+MM
      IF (II.GT.NN) II=II-NN
      IF (IFAC.EQ.9) GO TO 900
      IJ=II+MM
      IF (IJ.GT.NN) IJ=IJ-NN
      IK=IJ+MM
      IF (IK.GT.NN) IK=IK-NN
      IL=IK+MM
      IF (IL.GT.NN) IL=IL-NN
      IM=IL+MM
      IF (IM.GT.NN) IM=IM-NN
      IN=IM+MM
      IF (IN.GT.NN) IN=IN-NN
      IO=IN+MM
      IF (IO.GT.NN) IO=IO-NN
      IP=IO+MM
      IF (IP.GT.NN) IP=IP-NN
      IF (IFAC.EQ.16) GO TO 1600
*
*     ERROR EXIT - ILLEGAL FACTOR
*     ----------
      IERR=2
      GO TO 2000
*
C
C     CODING FOR FACTOR 2
C     -------------------
  200 CONTINUE
      DO 220 L=1,M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 210 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      A(IB+I) = ( AA - AB ) * FACT
      A(IA+I) = ( AA + AB ) * FACT
      BA = B(IA+I)
      BB = B(IB+I)
      B(IB+I) = ( BA - BB ) * FACT
      B(IA+I) = ( BA + BB ) * FACT
      I = I + JUMP
  210 CONTINUE
      IX=IB+INC
      IB=IA+INC
      IA=IX
  220 CONTINUE
      GO TO 1000
C
C     CODING FOR FACTOR 3
C     -------------------
  300 CONTINUE
      Z3=SIN60
      IF (MU.EQ.2) Z3=-Z3
*
      DO 340 L=1,M
*     --------------- REAL PARTS ---------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 310 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      T1 = AB + AC
      A(IC+I) = Z3 * ( AB - AC )
      A(IB+I) = AA - 0.5 * T1
      A(IA+I) = ( AA + T1 ) * FACT
      I = I + JUMP
  310 CONTINUE
*     ------------- IMAGINARY PARTS ------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 320 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      T1 = BB + BC
      B(IC+I) = Z3 * ( BB - BC )
      B(IB+I) = BA - 0.5 * T1
      B(IA+I) = ( BA + T1 ) * FACT
      I = I + JUMP
  320 CONTINUE
*     --------------- FOLDING ----------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 330 J = 1 , LOT
      T2 = B(IC+I)
      T3 = B(IB+I)
      T4 = A(IC+I)
      T1 = A(IB+I)
      A(IB+I) = ( T1 - T2 ) * FACT
      A(IC+I) = ( T1 + T2 ) * FACT
      B(IC+I) = ( T3 - T4 ) * FACT
      B(IB+I) = ( T3 + T4 ) * FACT
      I = I + JUMP
  330 CONTINUE
      IX=IC+INC
      IC=IB+INC
      IB=IA+INC
      IA=IX
  340 CONTINUE
      GO TO 1000
C
C     CODING FOR FACTOR 4
C     -------------------
  400 CONTINUE
      Z4=1.0
      IF (MU.EQ.3) Z4=-Z4
*
      DO 440 L=1,M
*     --------------- REAL PARTS ---------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 410 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      T1 = AA + AC
      T2 = AB + AD
      A(ID+I) = Z4 * ( AB - AD )
      A(IB+I) = AA - AC
      A(IA+I) = ( T1 + T2 ) * FACT
      A(IC+I) = ( T1 - T2 ) * FACT
      I = I + JUMP
  410 CONTINUE
*     ------------- IMAGINARY PARTS ------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 420 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      BD = B(ID+I)
      T1 = BA + BC
      T2 = BB + BD
      B(ID+I) = Z4 * ( BB - BD )
      B(IB+I) = BA - BC
      B(IA+I) = ( T1 + T2 ) * FACT
      B(IC+I) = ( T1 - T2 ) * FACT
      I = I + JUMP
  420 CONTINUE
*     --------------- FOLDING ----------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 430 J = 1 , LOT
      T2 = B(ID+I)
      T3 = B(IB+I)
      T4 = A(ID+I)
      T1 = A(IB+I)
      A(IB+I) = ( T1 - T2 ) * FACT
      A(ID+I) = ( T1 + T2 ) * FACT
      B(ID+I) = ( T3 - T4 ) * FACT
      B(IB+I) = ( T3 + T4 ) * FACT
      I = I + JUMP
  430 CONTINUE
      IX=ID+INC
      ID=IC+INC
      IC=IB+INC
      IB=IA+INC
      IA=IX
  440 CONTINUE
      GO TO 1000
C
C     CODING FOR FACTOR 5
C     -------------------
  500 CONTINUE
      IF (MU.EQ.1) THEN
        Z1=QRT5
        Z2=SIN72
        Z3=SIN36
      ELSE IF (MU.EQ.2) THEN
        Z1=-QRT5
        Z2=SIN36
        Z3=-SIN72
      ELSE IF (MU.EQ.3) THEN
        Z1=-QRT5
        Z2=-SIN36
        Z3=SIN72
      ELSE IF (MU.EQ.4) THEN
        Z1=QRT5
        Z2=-SIN72
        Z3=-SIN36
      ENDIF
*
      DO 540 L=1,M
*     --------------- REAL PARTS ---------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 510 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      AE = A(IE+I)
      T1 = AB + AE
      T2 = AC + AD
      T3 = T1 + T2
      T4 = Z1 * ( T1 - T2 )
      T1 = AB - AE
      T2 = AC - AD
      A(ID+I) = Z3 * T1 - Z2 * T2
      A(IE+I) = Z2 * T1 + Z3 * T2
      T1 = AA - 0.25 * T3
      A(IA+I) = ( AA + T3 ) * FACT
      A(IB+I) = T1 + T4
      A(IC+I) = T1 - T4
      I = I + JUMP
  510 CONTINUE
*     ------------- IMAGINARY PARTS ------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 520 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      BD = B(ID+I)
      BE = B(IE+I)
      T1 = BB + BE
      T2 = BC + BD
      T3 = T1 + T2
      T4 = Z1 * ( T1 - T2 )
      T1 = BB - BE
      T2 = BC - BD
      B(ID+I) = Z3 * T1 - Z2 * T2
      B(IE+I) = Z2 * T1 + Z3 * T2
      T1 = BA - 0.25 * T3
      B(IA+I) = ( BA + T3 ) * FACT
      B(IB+I) = T1 + T4
      B(IC+I) = T1 - T4
      I = I + JUMP
  520 CONTINUE
*     --------------- FOLDING ----------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 530 J = 1 , LOT
      T2 = B(IE+I)
      T3 = B(IB+I)
      T4 = A(IE+I)
      T1 = A(IB+I)
      A(IB+I) = ( T1 - T2 ) * FACT
      A(IE+I) = ( T1 + T2 ) * FACT
      B(IE+I) = ( T3 - T4 ) * FACT
      B(IB+I) = ( T3 + T4 ) * FACT
      T2 = B(ID+I)
      T3 = B(IC+I)
      T4 = A(ID+I)
      T1 = A(IC+I)
      A(IC+I) = ( T1 - T2 ) * FACT
      A(ID+I) = ( T1 + T2 ) * FACT
      B(ID+I) = ( T3 - T4 ) * FACT
      B(IC+I) = ( T3 + T4 ) * FACT
      I = I + JUMP
  530 CONTINUE
      IX=IE+INC
      IE=ID+INC
      ID=IC+INC
      IC=IB+INC
      IB=IA+INC
      IA=IX
  540 CONTINUE
      GO TO 1000
C
C     CODING FOR FACTOR 7
C     -------------------
  700 CONTINUE
      TWOPI=4.0*ASIN(1.0)
      ANGLE=FLOAT(MU)*TWOPI/7.0
      Z1=COS(ANGLE) 
      Z2=COS(2.0*ANGLE)
      Z3=COS(3.0*ANGLE)
      Z4=SIN(ANGLE) 
      Z5=SIN(2.0*ANGLE)
      Z6=SIN(3.0*ANGLE)
*
      DO 740 L=1,M
*     --------------- REAL PARTS ---------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 710 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      AE = A(IE+I)
      AF = A(IF+I)
      AG = A(IG+I)
      T1 = AB + AG
      AG = AB - AG
      T2 = AC + AF
      AF = AC - AF
      T3 = AD + AE
      AE = AD - AE
      T4 = AA - 0.5 * T3 
      A(IA+I) = ( AA + T1 + T2 + T3 ) * FACT
      T1 = T1 - T3
      T2 = T2 - T3
      A(IB+I) = T4 + Z1 * T1 + Z2 * T2
      A(IC+I) = T4 + Z2 * T1 + Z3 * T2
      A(ID+I) = T4 + Z3 * T1 + Z1 * T2
      A(IE+I) = Z6 * AG - Z4 * AF + Z5 * AE
      A(IF+I) = Z5 * AG - Z6 * AF - Z4 * AE
      A(IG+I) = Z4 * AG + Z5 * AF + Z6 * AE
      I = I + JUMP
  710 CONTINUE
*     ------------- IMAGINARY PARTS ------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 720 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      BD = B(ID+I)
      BE = B(IE+I)
      BF = B(IF+I)
      BG = B(IG+I)
      T1 = BB + BG
      BG = BB - BG
      T2 = BC + BF
      BF = BC - BF
      T3 = BD + BE
      BE = BD - BE
      T4 = BA - 0.5 * T3 
      B(IA+I) = ( BA + T1 + T2 + T3 ) * FACT
      T1 = T1 - T3
      T2 = T2 - T3
      B(IB+I) = T4 + Z1 * T1 + Z2 * T2
      B(IC+I) = T4 + Z2 * T1 + Z3 * T2
      B(ID+I) = T4 + Z3 * T1 + Z1 * T2
      B(IE+I) = Z6 * BG - Z4 * BF + Z5 * BE
      B(IF+I) = Z5 * BG - Z6 * BF - Z4 * BE
      B(IG+I) = Z4 * BG + Z5 * BF + Z6 * BE
      I = I + JUMP
  720 CONTINUE
*     --------------- FOLDING ----------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 730 J = 1 , LOT
      T2 = B(IG+I)
      T3 = B(IB+I)
      T4 = A(IG+I)
      T1 = A(IB+I)
      A(IB+I) = ( T1 - T2 ) * FACT
      A(IG+I) = ( T1 + T2 ) * FACT
      B(IG+I) = ( T3 - T4 ) * FACT
      B(IB+I) = ( T3 + T4 ) * FACT
      T2 = B(IF+I)
      T3 = B(IC+I)
      T4 = A(IF+I)
      T1 = A(IC+I)
      A(IC+I) = ( T1 - T2 ) * FACT
      A(IF+I) = ( T1 + T2 ) * FACT
      B(IF+I) = ( T3 - T4 ) * FACT
      B(IC+I) = ( T3 + T4 ) * FACT
      T2 = B(IE+I)
      T3 = B(ID+I)
      T4 = A(IE+I)
      T1 = A(ID+I)
      A(ID+I) = ( T1 - T2 ) * FACT
      A(IE+I) = ( T1 + T2 ) * FACT
      B(IE+I) = ( T3 - T4 ) * FACT
      B(ID+I) = ( T3 + T4 ) * FACT
      I = I + JUMP
  730 CONTINUE
      IX=IG+INC
      IG=IF+INC
      IF=IE+INC
      IE=ID+INC
      ID=IC+INC
      IC=IB+INC
      IB=IA+INC
      IA=IX
  740 CONTINUE
      GO TO 1000
C
C     CODING FOR FACTOR 8
C     -------------------
  800 CONTINUE
      Z1=1.0
      IF (MU.EQ.3.OR.MU.EQ.7) Z1=-Z1
      Z2=SIN45
      IF (MU.EQ.3.OR.MU.EQ.5) Z2=-Z2
      Z3=Z1*Z2
*
      DO 840 L=1,M
*     --------------- REAL PARTS ---------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 810 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      AE = A(IE+I)
      AF = A(IF+I)
      AG = A(IG+I)
      AH = A(IH+I)
      T1 = AA - AE
      AA = AA + AE
      T2 = AB - AF
      AE = AB + AF
      T3 = Z1 * ( AC - AG )
      AC = AC + AG
      T4 = AD -AH
      AG = AD +AH
      T5 = Z2 * ( T2 - T4 )
      T6 = Z3 * ( T2 + T4 )
      A(IB+I) = T1 + T5
      A(ID+I) = T1 - T5
      A(IF+I) = T6 - T3
      A(IH+I) = T6 + T3
      T1 = AA + AC
      A(IC+I) = AA - AC
      T2 = AE + AG
      A(IG+I) = Z1 * ( AE - AG )
      A(IA+I) = ( T1 + T2 ) * FACT
      A(IE+I) = ( T1 - T2 ) * FACT
      I = I + JUMP
  810 CONTINUE
*     ------------- IMAGINARY PARTS ------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 820 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      BD = B(ID+I)
      BE = B(IE+I)
      BF = B(IF+I)
      BG = B(IG+I)
      BH = B(IH+I)
      T1 = BA - BE
      BA = BA + BE
      T2 = BB - BF
      BE = BB + BF
      T3 = Z1 * ( BC - BG )
      BC = BC + BG
      T4 = BD - BH
      BG = BD + BH
      T5 = Z2 * ( T2 - T4 )
      T6 = Z3 * ( T2 + T4 )
      B(IB+I) = T1 + T5
      B(ID+I) = T1 - T5
      B(IF+I) = T6 - T3
      B(IH+I) = T6 + T3
      T1 = BA + BC
      B(IC+I) = BA - BC
      T2 = BE + BG
      B(IG+I) = Z1 * ( BE - BG )
      B(IA+I) = ( T1 + T2 ) * FACT
      B(IE+I) = ( T1 - T2 ) * FACT
      I = I + JUMP
  820 CONTINUE
*     --------------- FOLDING ----------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 830 J = 1 , LOT
      T2 = B(IH+I)
      T3 = B(IB+I)
      T4 = A(IH+I)
      T1 = A(IB+I)
      A(IB+I) = ( T1 - T2 ) * FACT
      A(IH+I) = ( T1 + T2 ) * FACT
      B(IH+I) = ( T3 - T4 ) * FACT
      B(IB+I) = ( T3 + T4 ) * FACT
      T2 = B(IG+I)
      T3 = B(IC+I)
      T4 = A(IG+I)
      T1 = A(IC+I)
      A(IC+I) = ( T1 - T2 ) * FACT
      A(IG+I) = ( T1 + T2 ) * FACT
      B(IG+I) = ( T3 - T4 ) * FACT
      B(IC+I) = ( T3 + T4 ) * FACT
      T2 = B(IF+I)
      T3 = B(ID+I)
      T4 = A(IF+I)
      T1 = A(ID+I)
      A(ID+I) = ( T1 - T2 ) * FACT
      A(IF+I) = ( T1 + T2 ) * FACT
      B(IF+I) = ( T3 - T4 ) * FACT
      B(ID+I) = ( T3 + T4 ) * FACT
      I = I + JUMP
  830 CONTINUE
      IX=IH+INC
      IH=IG+INC
      IG=IF+INC
      IF=IE+INC
      IE=ID+INC
      ID=IC+INC
      IC=IB+INC
      IB=IA+INC
      IA=IX
  840 CONTINUE
      GO TO 1000
C
C     CODING FOR FACTOR 9
C     -------------------
  900 CONTINUE
      Z1=SIN60
      IF (MOD(MU,3).EQ.2) Z1=-Z1
      TWOPI=4.0*ASIN(1.0)
      ANGLE=FLOAT(MU)*TWOPI/9.0
      Z2=COS(ANGLE) 
      Z3=SIN(ANGLE) 
      Z4=COS(2.0*ANGLE)
      Z5=SIN(2.0*ANGLE)
      Z6=Z1*Z2
      Z7=Z1*Z3
      Z8=Z1*Z4
      Z9=Z1*Z5
*
      DO 940 L=1,M
*     --------------- REAL PARTS ---------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 910 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      AE = A(IE+I)
      AF = A(IF+I)
      AG = A(IG+I)
      AH = A(IH+I)
      AI = A(II+I)
      T1 = AD + AG
      AD = Z1 * ( AD - AG )
      AG = AA - 0.5 * T1
      AA = AA + T1
      T2 = AE + AH
      AE = AE - AH
      AH = AB - 0.5 * T2
      AB = AB + T2
      T1 = AC + AF
      T2 = AC - AF
      T3 = AI - 0.5 * T1 
      AI = AI + T1
      T1 = AH + T3
      AH = AH - T3
      T3 = AE - T2
      AE = AE + T2
      T4 = Z2 * T1 - Z7 * T3
      T5 = Z4 * T1 + Z9 * T3
      A(IC+I) = AG + T5
      T1 = T4 + T5
      T2 = AG + T4
      T1 = AG - T1
      T3 = AB + AI
      A(IG+I) = Z1 * ( AB - AI )
      A(IB+I) = T2
      T4 = AA - 0.5 * T3 
      A(IA+I) = ( AA + T3 ) * FACT
      T2 = Z3 * AH + Z6 * AE
      T3 = Z5 * AH - Z8 * AE
      A(IE+I) = T1
      T5 = T2 - T3
      A(II+I) = AD + T2
      A(IF+I) = AD - T5
      A(IH+I) = T3 - AD
      A(ID+I) = T4
      I = I + JUMP
  910 CONTINUE
*     ------------- IMAGINARY PARTS ------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 920 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      BD = B(ID+I)
      BE = B(IE+I)
      BF = B(IF+I)
      BG = B(IG+I)
      BH = B(IH+I)
      BI = B(II+I)
      T1 = BD + BG
      BD = Z1 * ( BD - BG )
      BG = BA - 0.5 * T1
      BA = BA + T1
      T2 = BE + BH
      BE = BE - BH
      BH = BB - 0.5 * T2
      BB = BB + T2
      T1 = BC + BF
      T2 = BC - BF
      T3 = BI - 0.5 * T1 
      BI = BI + T1
      T1 = BH + T3
      BH = BH - T3
      T3 = BE - T2
      BE = BE + T2
      T4 = Z2 * T1 - Z7 * T3
      T5 = Z4 * T1 + Z9 * T3
      B(IC+I) = BG + T5
      T1 = T4 + T5
      T2 = BG + T4
      T1 = BG - T1
      T3 = BB + BI
      B(IG+I) = Z1 * ( BB - BI )
      B(IB+I) = T2
      T4 = BA - 0.5 * T3 
      B(IA+I) = ( BA + T3 ) * FACT
      T2 = Z3 * BH + Z6 * BE
      T3 = Z5 * BH - Z8 * BE
      B(IE+I) = T1
      T5 = T2 - T3
      B(II+I) = BD + T2
      B(IF+I) = BD - T5
      B(IH+I) = T3 - BD
      B(ID+I) = T4
      I = I + JUMP
  920 CONTINUE
*     --------------- FOLDING ----------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 930 J = 1 , LOT
      T2 = B(II+I)
      T3 = B(IB+I)
      T4 = A(II+I)
      T1 = A(IB+I)
      A(IB+I) = ( T1 - T2 ) * FACT
      A(II+I) = ( T1 + T2 ) * FACT
      B(II+I) = ( T3 - T4 ) * FACT
      B(IB+I) = ( T3 + T4 ) * FACT
      T2 = B(IH+I)
      T3 = B(IC+I)
      T4 = A(IH+I)
      T1 = A(IC+I)
      A(IC+I) = ( T1 - T2 ) * FACT
      A(IH+I) = ( T1 + T2 ) * FACT
      B(IH+I) = ( T3 - T4 ) * FACT
      B(IC+I) = ( T3 + T4 ) * FACT
      T2 = B(IG+I)
      T3 = B(ID+I)
      T4 = A(IG+I)
      T1 = A(ID+I)
      A(ID+I) = ( T1 - T2 ) * FACT
      A(IG+I) = ( T1 + T2 ) * FACT
      B(IG+I) = ( T3 - T4 ) * FACT
      B(ID+I) = ( T3 + T4 ) * FACT
      T2 = B(IF+I)
      T3 = B(IE+I)
      T4 = A(IF+I)
      T1 = A(IE+I)
      A(IE+I) = ( T1 - T2 ) * FACT
      A(IF+I) = ( T1 + T2 ) * FACT
      B(IF+I) = ( T3 - T4 ) * FACT
      B(IE+I) = ( T3 + T4 ) * FACT
      I = I + JUMP
  930 CONTINUE
      IX=II+INC
      II=IH+INC
      IH=IG+INC
      IG=IF+INC
      IF=IE+INC
      IE=ID+INC
      ID=IC+INC
      IC=IB+INC
      IB=IA+INC
      IA=IX
  940 CONTINUE
      GO TO 1000
C
C     CODING FOR FACTOR 16
C     --------------------
 1600 CONTINUE
      Z1=1.0
      IF (MOD(MU,4).EQ.3) Z1=-Z1
      ANGLE = FLOAT(MU) * 0.25 * ASIN(1.0)
      Z2=COS(ANGLE) 
      Z3=SIN(ANGLE) 
      Z4=SIN45
      MMU=MOD(MU,8) 
      IF (MMU.EQ.3.OR.MMU.EQ.5) Z4=-Z4
      Z5=Z1*Z4
      Z6=Z1*Z3
      Z7=Z1*Z2
*
      DO 1660 L=1,M 
*     --------------- REAL PARTS ---------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1610 J = 1 , LOT

      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      AE = A(IE+I)
      AF = A(IF+I)
      AG = A(IG+I)
      AH = A(IH+I)
      AI = A(II+I)
      AJ = A(IJ+I)
      AK = A(IK+I)
      AL = A(IL+I)
      AM = A(IM+I)
      AN = A(IN+I)
      AO = A(IO+I)
      AP = A(IP+I)

      T1 = AA + AI
      T2 = AE + AM
      AE = Z1 * ( AE - AM )
      AM = AA - AI
      AA = T1 + T2
      AI = T1 - T2
      T2 = AF + AN
      AF = AF - AN
      T1 = AB + AJ
      AN = AB - AJ
      AB = T1 + T2
      AJ = T1 - T2
      T1 = AC + AK
      T2 = AC - AK
      T3 = AG + AO
      T4 = AG - AO
      AC = T1 + T3
      AG = Z1 * ( T1 - T3 )
      AK = Z4 * ( T2 - T4 )
      AO = Z5 * ( T2 + T4 )
      T1 = AD + AL
      T2 = AD - AL
      T3 = AH + AP
      T4 = AH - AP
      AD = T1 + T3
      AL = T1 - T3
      T1 = AN - T4
      T3 = AN + T4
      T4 = AF + T2
      T2 = AF - T2
      AH = Z2 * T1 - Z6 * T2
      AF = Z6 * T1 + Z2 * T2
      AN = Z7 * T3 - Z3 * T4
      AP = Z3 * T3 + Z7 * T4

      A(IE+I) = AA - AC
      T1 = AA + AC
      T2 = AB + AD
      A(IM+I) = Z1 * ( AB - AD )
      A(IA+I) = ( T1 + T2 ) * FACT
      A(II+I) = ( T1 - T2 ) * FACT
      T1 = AM + AK
      T2 = AM - AK
      A(IB+I) = T1 + AH
      A(IH+I) = T1 - AH
      A(ID+I) = T2 + AF
      A(IF+I) = T2 - AF
      T1 = AE + AO
      T2 = AE - AO
      A(IJ+I) = AP - T1
      A(IP+I) = AP + T1
      A(IL+I) = AN + T2
      A(IN+I) = AN - T2
      T1 = Z4 * ( AJ - AL )
      T2 = Z5 * ( AJ + AL )
      A(IK+I) = T2 - AG
      A(IO+I) = T2 + AG
      A(IC+I) = AI + T1
      A(IG+I) = AI - T1

      I = I + JUMP

 1610 CONTINUE
*     ------------- IMAGINARY PARTS ------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1630 J = 1 , LOT

      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      BD = B(ID+I)
      BE = B(IE+I)
      BF = B(IF+I)
      BG = B(IG+I)
      BH = B(IH+I)
      BI = B(II+I)
      BJ = B(IJ+I)
      BK = B(IK+I)
      BL = B(IL+I)
      BM = B(IM+I)
      BN = B(IN+I)
      BO = B(IO+I)
      BP = B(IP+I)

      T1 = BA + BI
      T2 = BE + BM
      BE = Z1 * ( BE - BM )
      BM = BA - BI
      BA = T1 + T2
      BI = T1 - T2
      T2 = BF + BN
      BF = BF - BN
      T1 = BB + BJ
      BN = BB - BJ
      BB = T1 + T2
      BJ = T1 - T2
      T1 = BC + BK
      T2 = BC - BK
      T3 = BG + BO
      T4 = BG - BO
      BC = T1 + T3
      BG = Z1 * ( T1 - T3 )
      BK = Z4 * ( T2 - T4 )
      BO = Z5 * ( T2 + T4 )
      T1 = BD + BL
      T2 = BD - BL
      T3 = BH + BP
      T4 = BH - BP
      BD = T1 + T3
      BL = T1 - T3
      T1 = BN - T4
      T3 = BN + T4
      T4 = BF + T2
      T2 = BF - T2
      BH = Z2 * T1 - Z6 * T2
      BF = Z6 * T1 + Z2 * T2
      BN = Z7 * T3 - Z3 * T4
      BP = Z3 * T3 + Z7 * T4

      B(IE+I) = BA - BC
      T1 = BA + BC
      T2 = BB + BD
      B(IM+I) = Z1 * ( BB - BD )
      B(IA+I) = ( T1 + T2 ) * FACT
      B(II+I) = ( T1 - T2 ) * FACT
      T1 = BM + BK
      T2 = BM - BK
      B(IB+I) = T1 + BH
      B(IH+I) = T1 - BH
      B(ID+I) = T2 + BF
      B(IF+I) = T2 - BF
      T1 = BE + BO
      T2 = BE - BO
      B(IJ+I) = BP - T1
      B(IP+I) = BP + T1
      B(IL+I) = BN + T2
      B(IN+I) = BN - T2
      T1 = Z4 * ( BJ - BL )
      T2 = Z5 * ( BJ + BL )
      B(IK+I) = T2 - BG
      B(IO+I) = T2 + BG
      B(IC+I) = BI + T1
      B(IG+I) = BI - T1

      I = I + JUMP

 1630 CONTINUE
*     --------------- FOLDING ----------------
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1650 J = 1 , LOT
      T2 = B(IP+I)
      T3 = B(IB+I)
      T4 = A(IP+I)
      T1 = A(IB+I)
      A(IB+I) = ( T1 - T2 ) * FACT
      A(IP+I) = ( T1 + T2 ) * FACT
      B(IP+I) = ( T3 - T4 ) * FACT
      B(IB+I) = ( T3 + T4 ) * FACT
      T2 = B(IO+I)
      T3 = B(IC+I)
      T4 = A(IO+I)
      T1 = A(IC+I)
      A(IC+I) = ( T1 - T2 ) * FACT
      A(IO+I) = ( T1 + T2 ) * FACT
      B(IO+I) = ( T3 - T4 ) * FACT
      B(IC+I) = ( T3 + T4 ) * FACT
      T2 = B(IN+I)
      T3 = B(ID+I)
      T4 = A(IN+I)
      T1 = A(ID+I)
      A(ID+I) = ( T1 - T2 ) * FACT
      A(IN+I) = ( T1 + T2 ) * FACT
      B(IN+I) = ( T3 - T4 ) * FACT
      B(ID+I) = ( T3 + T4 ) * FACT
      T2 = B(IM+I)
      T3 = B(IE+I)
      T4 = A(IM+I)
      T1 = A(IE+I)
      A(IE+I) = ( T1 - T2 ) * FACT
      A(IM+I) = ( T1 + T2 ) * FACT
      B(IM+I) = ( T3 - T4 ) * FACT
      B(IE+I) = ( T3 + T4 ) * FACT
      T2 = B(IL+I)
      T3 = B(IF+I)
      T4 = A(IL+I)
      T1 = A(IF+I)
      A(IF+I) = ( T1 - T2 ) * FACT
      A(IL+I) = ( T1 + T2 ) * FACT
      B(IL+I) = ( T3 - T4 ) * FACT
      B(IF+I) = ( T3 + T4 ) * FACT
      T2 = B(IK+I)
      T3 = B(IG+I)
      T4 = A(IK+I)
      T1 = A(IG+I)
      A(IG+I) = ( T1 - T2 ) * FACT
      A(IK+I) = ( T1 + T2 ) * FACT
      B(IK+I) = ( T3 - T4 ) * FACT
      B(IG+I) = ( T3 + T4 ) * FACT
      T2 = B(IJ+I)
      T3 = B(IH+I)
      T4 = A(IJ+I)
      T1 = A(IH+I)
      A(IH+I) = ( T1 - T2 ) * FACT
      A(IJ+I) = ( T1 + T2 ) * FACT
      B(IJ+I) = ( T3 - T4 ) * FACT
      B(IH+I) = ( T3 + T4 ) * FACT
      I = I + JUMP
 1650 CONTINUE
      IX=IP+INC
      IP=IO+INC
      IO=IN+INC
      IN=IM+INC
      IM=IL+INC
      IL=IK+INC
      IK=IJ+INC
      IJ=II+INC
      II=IH+INC
      IH=IG+INC
      IG=IF+INC
      IF=IE+INC
      IE=ID+INC
      ID=IC+INC
      IC=IB+INC
      IB=IA+INC
      IA=IX
 1660 CONTINUE
C
 1000 CONTINUE
*
      IERR=0
 2000 CONTINUE
      RETURN
      END 
*DECK FORPFA
*        SUBROUTINE 'RPFA2'
*        SELF-SORTING IN-PLACE PRIME FACTOR FFT
*         (REAL TO HALF-COMPLEX)
*
*        CALL RPFA2(A,B,IFAX,ILIST,INC1,JUMP1,INC2,JUMP2,
*       +           N,LOT,IFORM,IERR)
*
*        A IS THE INPUT ARRAY (GRIDPOINT VALUES)
*        B IS THE OUTPUT ARRAY (SPECTRAL COEFFICIENTS)
*        IFAX IS A LIST OF FACTORS OF N 
*              IFAX(1) CONTAINS THE NUMBER OF FACTORS (NFAX)
*              THE FACTORS FOLLOW IN IFAX(2) TO IFAX(NFAX+1)
*        ILIST IS A LIST OF ADDRESSES USED IN THE 'FOLDING' 
*          STAGE OF THE ALGORITHM
*        INC1 IS THE INCREMENT WITHIN THE INPUT DATA VECTORS
*        JUMP1 IS THE INCREMENT BETWEEN THE INPUT DATA VECTORS
*        INC2 IS THE INCREMENT WITHIN THE OUTPUT DATA VECTORS
*        JUMP2 IS THE INCREMENT BETWEEN THE OUTPUT DATA VECTORS
*        N IS THE LENGTH OF THE TRANSFORMS
*        LOT IS THE NUMBER OF TRANSFORMS
*        IFORM = +1 FOR CONVENTIONAL ORDERING OF SPECTRAL COEFFS
*                -1 FOR UNCONVENTIONAL (FAST-IN-PLACE) ORDERING
*        IERR IS AN ERROR INDICATOR: 
*                  0 - TRANSFORM COMPLETED WITHOUT ERROR
*                  1 - FACTORS ARE NOT MUTUALLY PRIME
*                  2 - ILLEGAL FACTOR
*
*        FORTRAN VERSION WRITTEN BY CLIVE TEMPERTON
*        RECHERCHE EN PREVISION NUMERIQUE -- OCTOBER 1986
*
*        MODIFIED BY CLIVE TEMPERTON MAY 1987 -- TO GUARANTEE
*        CORRECT IN-PLACE WORKING OF FORTRAN VERSION
*
*        MODIFIED BY BERNARD DUGAS DECEMBER 1992 -- USE "IMPLICIT NONE"
*        AND USE VECTOR TEMPORARIES WITHIN LOOPS
*
*----------------------------------------------------------------------
*
*        FOR A MATHEMATICAL DEVELOPMENT OF THE ALGORITHM USED,
*        SEE: 
*
*        C TEMPERTON : "IMPLEMENTATION OF A SELF-SORTING IN-PLACE
*          PRIME FACTOR FFT ALGORITHM", 
*          JOURNAL OF COMPUTATIONAL PHYSICS VOL 58 (1985), 283-299.
*
*        C TEMPERTON : "A SELF-SORTING IN-PLACE PRIME FACTOR
*          REAL/HALF-COMPLEX FFT ALGORITHM",
*          TO APPEAR IN JOURNAL OF COMPUTATIONAL PHYSICS
*
*        C TEMPERTON : "A NEW SET OF MINIMUM-ADD SMALL-N DFT MODULES",
*          TO APPEAR IN JOURNAL OF COMPUTATIONAL PHYSICS
*
*----------------------------------------------------------------------
*
      SUBROUTINE RPFA2(A,B,IFAX,ILIST,INC1,JUMP1,INC2,JUMP2, 
     *                 N,LOT,IFORM,IERR) 
*
      IMPLICIT    none
*
      INTEGER     INC1,JUMP1,INC2,JUMP2,N,LOT,IFORM,IERR
      REAL        A(*),B(*)
      INTEGER     IFAX(*),ILIST(*)
*
      INTEGER     MAXLOT
      PARAMETER ( MAXLOT = lot_maximum )
*
      REAL       V1(MAXLOT+1), V2(MAXLOT+1), V3(MAXLOT+1), V4(MAXLOT+1), 
     *           V5(MAXLOT+1), V6(MAXLOT+1), V7(MAXLOT+1), V8(MAXLOT+1), 
     *           V9(MAXLOT+1),V10(MAXLOT+1),V11(MAXLOT+1),V12(MAXLOT+1),
     *          V13(MAXLOT+1),V14(MAXLOT+1),V15(MAXLOT+1),V16(MAXLOT+1) 
*
      REAL        AA,AB,AC,AD,AE,AF,AG,AH,AI,AJ,AK,AL,AM,AN,AO,AP
      REAL        S,S1,S2,S3,S4,S5,S6,S7,S8,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9
      REAL        T1,T2,T3,T4,T5,T6
*
      INTEGER     NBLOX,I,J,K,L,LV,LV1,M,MM,MU,NN,NB
      INTEGER     IA,IB,IC,ID,IE,IF,IG,IH,II,IJ,IK,IL,IM,IN,IO,IP,IQ
      INTEGER     JA,JB,JC,JD,JE,JF,JG,JH,JI,JJ,JK,JL,JM,JN,JO,JP
      INTEGER     NFAX,IFAC,ISHIFT,IX

      REAL        SIN60,QRT5, SIN36,SIN72,SIN45
      REAL        COS26,SIN26,COS51,SIN51
      REAL        COS77,SIN77,COS20,SIN20
      REAL        COS40,SIN40,COS80,SIN80
      REAL        COS23,SIN23
*
      DATA   SIN60/0.866025403784437/
      DATA   QRT5 /0.559016994374947/
      DATA   SIN36/0.587785252292473/,SIN72/0.951056516295154/
      DATA   COS26/0.900968867902419/,SIN26/0.433883739117558/
      DATA   COS51/0.623489801858734/,SIN51/0.781831482468030/
      DATA   COS77/0.222520933956314/,SIN77/0.974927912181824/
      DATA   SIN45/0.707106781186548/
      DATA   SIN20/0.342020143325669/,COS20/0.939692620785908/,
     *       SIN40/0.642787609686539/,COS40/0.766044443118978/,
     *       SIN80/0.984807753012208/,COS80/0.173648177666930/
      DATA   SIN23/0.382683432365090/,COS23/0.923879532511287/
C
      NFAX = IFAX(1)
C
      DO 1000 K = 1 , NFAX
C
      IFAC=IFAX(K+1)
      M=N/IFAC
      DO 100 J = 1 , IFAC
      MU=J
      MM=J*M
      IF (MOD(MM,IFAC).EQ.1) GO TO 110
  100 CONTINUE
*
*     ERROR EXIT (FACTORS NOT MUTUALLY PRIME)
*     ----------
      IERR = 1
      GO TO 2000
*
  110 CONTINUE
      MU = IFAC - MU
      MM = MM * INC1
      NN = N * INC1 
C
      IA=1
      IB=IA+MM
      IF (IFAC.EQ.2) GO TO 200
      IC=IB+MM
      IF (IC.GT.NN) IC=IC-NN
      IF (IFAC.EQ.3) GO TO 300
      ID=IC+MM
      IF (ID.GT.NN) ID=ID-NN
      IF (IFAC.EQ.4) GO TO 400
      IE=ID+MM
      IF (IE.GT.NN) IE=IE-NN
      IF (IFAC.EQ.5) GO TO 500
      IF=IE+MM
      IF (IF.GT.NN) IF=IF-NN
      IG=IF+MM
      IF (IG.GT.NN) IG=IG-NN
      IF (IFAC.EQ.7) GO TO 700
      IH=IG+MM
      IF (IH.GT.NN) IH=IH-NN
      IF (IFAC.EQ.8) GO TO 800
      II=IH+MM
      IF (II.GT.NN) II=II-NN
      IF (IFAC.EQ.9) GO TO 900
      IJ=II+MM
      IF (IJ.GT.NN) IJ=IJ-NN
      IK=IJ+MM
      IF (IK.GT.NN) IK=IK-NN
      IL=IK+MM
      IF (IL.GT.NN) IL=IL-NN
      IM=IL+MM
      IF (IM.GT.NN) IM=IM-NN
      IN=IM+MM
      IF (IN.GT.NN) IN=IN-NN
      IO=IN+MM
      IF (IO.GT.NN) IO=IO-NN
      IP=IO+MM
      IF (IP.GT.NN) IP=IP-NN
      IF (IFAC.EQ.16) GO TO 1600
*
*     ERROR EXIT (ILLEGAL FACTOR)
*     ----------
      IERR = 2
      GO TO 2000
C
C     CODING FOR FACTOR 2
C     -------------------
  200 CONTINUE
C
      DO 220 L = 1 , M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 210 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      A(IB+I) = AA - AB
      A(IA+I) = AA + AB
      I = I + JUMP1 
  210 CONTINUE
      IX=IB+INC1
      IB=IA+INC1
      IA=IX
  220 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 3
C     -------------------
  300 CONTINUE
      Z3=SIN60
      IF (MU.EQ.2) Z3=-Z3
C
      DO 320 L = 1 , M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 310 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      T1 = AB + AC
      A(IC+I) = Z3 * ( AB - AC )
      A(IB+I) = AA - .5 * T1
      A(IA+I) = AA + T1
      I = I + JUMP1 
  310 CONTINUE
      IX=IC+INC1
      IC=IB+INC1
      IB=IA+INC1
      IA=IX
  320 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 4
C     -------------------
  400 CONTINUE
      Z4=1.0
      IF (MU.EQ.3) Z4=-Z4
C
      DO 420 L = 1 , M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 410 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      T2 = AB + AD
      A(ID+I) = Z4 * ( AB - AD )
      A(IB+I) = AA - AC
      T1 = AA + AC
      A(IA+I) = T1 + T2
      A(IC+I) = T1 - T2
      I = I + JUMP1 
  410 CONTINUE
      IX=ID+INC1
      ID=IC+INC1
      IC=IB+INC1
      IB=IA+INC1
      IA=IX
  420 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 5
C     -------------------
  500 CONTINUE
      IF (MU.EQ.1) THEN
        Z1=QRT5
        Z2=SIN72
        Z3=SIN36
      ELSE IF (MU.EQ.2) THEN
        Z1=-QRT5
        Z2=SIN36
        Z3=-SIN72
      ELSE IF (MU.EQ.3) THEN
        Z1=-QRT5
        Z2=-SIN36
        Z3=SIN72
      ELSE IF (MU.EQ.4) THEN
        Z1=QRT5
        Z2=-SIN72
        Z3=-SIN36
      ENDIF
      DO 520 L=1,M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 510 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      AE = A(IE+I)
      T1 = AB + AE
      T2 = AC + AD
      T3 = AA - .25 * ( T1 + T2 ) 
      A(IA+I) =  AA + ( T1 + T2 )
      T4 = Z1 * ( T1 - T2 )
      T1 = AB - AE
      T2 = AC - AD
      A(IB+I) = T3 + T4
      A(IC+I) = T3 - T4
      A(ID+I) = Z3 * T1 - Z2 * T2
      A(IE+I) = Z2 * T1 + Z3 * T2
      I = I + JUMP1 
  510 CONTINUE
      IX = IE + INC1
      IE = ID + INC1
      ID = IC + INC1
      IC = IB + INC1
      IB = IA + INC1
      IA = IX
  520 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 7
C     -------------------
  700 CONTINUE
      IF (MU.EQ.1) THEN
        Z1=COS51
        Z2=-COS77
        Z3=-COS26
        Z4=SIN51
        Z5=SIN77
        Z6=SIN26
      ELSE IF (MU.EQ.2) THEN
        Z1=-COS77
        Z2=-COS26
        Z3=COS51
        Z4=SIN77
        Z5=-SIN26
        Z6=-SIN51
      ELSE IF (MU.EQ.3) THEN
        Z1=-COS26
        Z2=COS51
        Z3=-COS77
        Z4=SIN26
        Z5=-SIN51
        Z6=SIN77
      ELSE IF (MU.EQ.4) THEN
        Z1=-COS26
        Z2=COS51
        Z3=-COS77
        Z4=-SIN26
        Z5=SIN51
        Z6=-SIN77
      ELSE IF (MU.EQ.5) THEN
        Z1=-COS77
        Z2=-COS26
        Z3=COS51
        Z4=-SIN77
        Z5=SIN26
        Z6=SIN51
      ELSE IF (MU.EQ.6) THEN
        Z1=COS51
        Z2=-COS77
        Z3=-COS26
        Z4=-SIN51
        Z5=-SIN77
        Z6=-SIN26
      ENDIF
      DO 720 L=1,M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 710 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      AE = A(IE+I)
      AF = A(IF+I)
      AG = A(IG+I)
      T1 = AB + AG
      T2 = AC + AF
      T3 = AD + AE
      T4 = AB - AG
      T5 = AC - AF
      T6 = AD - AE
      A(IB+I) = AA + Z1 * T1 + Z2 * T2 + Z3 * T3
      A(IC+I) = AA + Z2 * T1 + Z3 * T2 + Z1 * T3
      A(ID+I) = AA + Z3 * T1 + Z1 * T2 + Z2 * T3
      A(IA+I) = AA +      T1 +      T2 +      T3
      A(IE+I) =      Z6 * T4 - Z4 * T5 + Z5 * T6
      A(IF+I) =      Z5 * T4 - Z6 * T5 - Z4 * T6
      A(IG+I) =      Z4 * T4 + Z5 * T5 + Z6 * T6
      I = I + JUMP1 
  710 CONTINUE
      IX = IG + INC1
      IG = IF + INC1
      IF = IE + INC1
      IE = ID + INC1
      ID = IC + INC1
      IC = IB + INC1
      IB = IA + INC1
      IA = IX
  720 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 8
C     -------------------
  800 CONTINUE
      IF (MU.EQ.1) THEN
        Z1=1.0
        Z2=SIN45
        Z3=SIN45
      ELSE IF (MU.EQ.3) THEN
        Z1=-1.0
        Z2=-SIN45
        Z3=SIN45
      ELSE IF (MU.EQ.5) THEN
        Z1=1.0
        Z2=-SIN45
        Z3=-SIN45
      ELSE IF (MU.EQ.7) THEN
        Z1=-1.0
        Z2=SIN45
        Z3=-SIN45
      ENDIF
      DO 820 L=1,M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 810 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      AE = A(IE+I)
      AF = A(IF+I)
      AG = A(IG+I)
      AH = A(IH+I)
      T1 = AA + AE
      T2 = AC + AG
      T3 = Z1 * ( AC - AG )
      T4 = T1 + T2
      A(IC+I) = T1 - T2
      T1 = AB - AF
      T2 = AD - AH
      T5 = Z2 * ( T1 - T2 )
      T6 = Z3 * ( T1 + T2 )
      T1 = AB + AF
      T2 = AD + AH
      A(IF+I) = T6 - T3
      A(IH+I) = T6 + T3
      T3 = AA - AE
      A(IB+I) = T3 + T5
      A(ID+I) = T3 - T5
      T6 = T1 + T2
      A(IG+I) = Z1 * ( T1 - T2 )
      A(IA+I) = T4 + T6
      A(IE+I) = T4 - T6
      I = I + JUMP1 
  810 CONTINUE
      IX = IH + INC1
      IH = IG + INC1
      IG = IF + INC1
      IF = IE + INC1
      IE = ID + INC1
      ID = IC + INC1
      IC = IB + INC1
      IB = IA + INC1
      IA = IX
  820 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 9
C     -------------------
  900 CONTINUE
      IF (MU.EQ.1) THEN
        Z1=SIN60
        Z2=COS40
        Z3=SIN40
        Z4=COS80
        Z5=SIN80
      ELSE IF (MU.EQ.2) THEN
        Z1=-SIN60
        Z2=COS80
        Z3=SIN80
        Z4=-COS20
        Z5=SIN20
      ELSE IF (MU.EQ.4) THEN
        Z1=SIN60
        Z2=-COS20
        Z3=SIN20
        Z4=COS40
        Z5=-SIN40
      ELSE IF (MU.EQ.5) THEN
        Z1=-SIN60
        Z2=-COS20
        Z3=-SIN20
        Z4=COS40
        Z5=SIN40
      ELSE IF (MU.EQ.7) THEN
        Z1=SIN60
        Z2=COS80
        Z3=-SIN80
        Z4=-COS20
        Z5=-SIN20
      ELSE IF (MU.EQ.8) THEN
        Z1=-SIN60
        Z2=COS40
        Z3=-SIN40
        Z4=COS80
        Z5=-SIN80
      ENDIF
      Z6=Z1*Z2
      Z7=Z1*Z3
      Z8=Z1*Z4
      Z9=Z1*Z5
      DO 920 L = 1 , M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 910 J = 1 , LOT
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      AE = A(IE+I)
      AF = A(IF+I)
      AG = A(IG+I)
      AH = A(IH+I)
      AI = A(II+I)
      T1 = AD + AG
      T2 = AA - .5 * T1 
      AD = Z1 * ( AD - AG )
      AA = AA + T1
      T1 = AE + AH
      T3 = AB - .5 * T1 
      T4 = AE - AH
      AG = AB + T1
      AE = Z2 * T3 - Z7 * T4
      AH = Z3 * T3 + Z6 * T4
      T1 = AF + AI
      T3 = AC - .5 * T1 
      T4 = AF - AI
      AC = AC + T1
      T1 = Z4 * T3 - Z9 * T4
      AI = Z5 * T3 + Z8 * T4
      T3 = AE + T1
      T4 = Z1 * ( AE - T1 )
      A(IB+I) = T2 + T3
      T2 = T2 - .5 * T3
      T1 = AH + AI
      AE = Z1 * ( AH - AI )
      A(II+I) = AD + T1
      T3 = AD - .5 * T1 
      T1 = AG + AC
      A(IG+I) = Z1 * ( AG - AC )
      A(ID+I) = AA - .5 * T1
      A(IA+I) = AA + T1
      A(IF+I) = T3 + T4
      A(IH+I) = T4 - T3
      A(IC+I) = T2 + AE
      A(IE+I) = T2 - AE
      I = I + JUMP1 
  910 CONTINUE
      IX = II + INC1
      II = IH + INC1
      IH = IG + INC1
      IG = IF + INC1
      IF = IE + INC1
      IE = ID + INC1
      ID = IC + INC1
      IC = IB + INC1
      IB = IA + INC1
      IA = IX
  920 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 16
C     --------------------
 1600 CONTINUE
      IF (MU.EQ.1) THEN
        Z1=1.0
        Z2=COS23
        Z3=SIN23
        Z4=SIN45
      ELSE IF (MU.EQ.3) THEN
        Z1=-1.0
        Z2=SIN23
        Z3=COS23
        Z4=-SIN45
      ELSE IF (MU.EQ.5) THEN
        Z1=1.0
        Z2=-SIN23
        Z3=COS23
        Z4=-SIN45
      ELSE IF (MU.EQ.7) THEN
        Z1=-1.0
        Z2=-COS23
        Z3=SIN23
        Z4=SIN45
      ELSE IF (MU.EQ.9) THEN
        Z1=1.0
        Z2=-COS23
        Z3=-SIN23
        Z4=SIN45
      ELSE IF (MU.EQ.11) THEN 
        Z1=-1.0
        Z2=-SIN23
        Z3=-COS23
        Z4=-SIN45
      ELSE IF (MU.EQ.13) THEN 
        Z1=1.0
        Z2=SIN23
        Z3=-COS23
        Z4=-SIN45
      ELSE IF (MU.EQ.15) THEN 
        Z1=-1.0
        Z2=COS23
        Z3=-SIN23
        Z4=SIN45
      ENDIF
      Z5=Z1*Z4
      Z6=Z1*Z3
      Z7=Z1*Z2
      DO 1620 L = 1 , M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1610 J = 1 , LOT
      AA = A(IA+I)
      AE = A(IE+I)
      AI = A(II+I)
      AM = A(IM+I)
      T1 = AA + AI
      T2 = AE + AM
      AI = AA - AI
      T3 = Z1 * ( AE - AM )
      AA = T1 + T2
      T4 = T1 - T2
      AB = A(IB+I)
      AF = A(IF+I)
      AJ = A(IJ+I)
      AN = A(IN+I)
      T1 = AB + AJ
      T2 = AF + AN
      AJ = AB - AJ
      AN = AF - AN
      AB = T1 + T2
      AF = T1 - T2
      AC = A(IC+I)
      AG = A(IG+I)
      AK = A(IK+I)
      AO = A(IO+I)
      T1 = AC - AK
      AK = AC + AK
      T2 = AG + AO
      AG = AG - AO
      AC = AK + T2
      AO = Z1 * ( AK - T2 )
      AK = Z4 * ( T1 - AG )
      AG = Z5 * ( T1 + AG )
      AD = A(ID+I)
      AH = A(IH+I)
      AL = A(IL+I)
      AP = A(IP+I)
      T1 = AD + AL
      AL = AD - AL
      T2 = AH + AP
      AP = AH - AP
      AD = T1 + T2
      AH = T1 - T2
      T1 = AJ + AP
      AP = AJ - AP
      T2 = AN + AL
      AL = AN - AL
      AJ = Z7 * T1 - Z3 * T2
      AN = Z3 * T1 + Z7 * T2
      T1 = AA + AC
      A(IE+I) = AA - AC
      T2 = AB + AD
      A(IM+I) = Z1 * ( AB - AD )
      T5 = AI + AK
      T6 = AI - AK
      A(IA+I) = T1 + T2
      A(II+I) = T1 - T2
      T1 = Z4 * ( AF - AH )
      T2 = Z5 * ( AF + AH )
      A(IC+I) = T4 + T1
      T1 = T4 - T1
      T4 = T2 + AO
      A(IK+I) = T2 - AO
      A(IO+I) = T4
      T2 = T3 + AG
      T4 = T3 - AG
      A(IG+I) = T1
      T1 = Z2 * AP - Z6 * AL
      T3 = Z6 * AP + Z2 * AL
      A(IB+I) = T5 + T1
      A(IH+I) = T5 - T1
      A(ID+I) = T6 + T3
      A(IF+I) = T6 - T3
      A(IL+I) = AJ + T4
      T1 = AJ - T4
      A(IJ+I) = AN - T2
      A(IP+I) = AN + T2
      A(IN+I) = T1
      I = I + JUMP1 
 1610 CONTINUE
      IX = IP + INC1
      IP = IO + INC1
      IO = IN + INC1
      IN = IM + INC1
      IM = IL + INC1
      IL = IK + INC1
      IK = IJ + INC1
      IJ = II + INC1
      II = IH + INC1
      IH = IG + INC1
      IG = IF + INC1
      IF = IE + INC1
      IE = ID + INC1
      ID = IC + INC1
      IC = IB + INC1
      IB = IA + INC1
      IA = IX
 1620 CONTINUE
C
 1000 CONTINUE
*
      NBLOX = 1 + (LOT-1)/MAXLOT
      LV1 = LOT - (NBLOX-1)*MAXLOT
*
*     K=0 AND K=N/2 
*     ------------- 
      JJ = INC2 + 1 
      S = 1.0 / FLOAT(N)
      IF (IFORM.EQ.+1) THEN
      DO 1005 L = 1 , LOT
      B(JJ) = 0.0
      JJ = JJ + JUMP2
 1005 CONTINUE
      ENDIF
      I = 1
      J = 1
      DO 1010 L = 1 , LOT
      B(J) = S * A(I)
      I = I + JUMP1 
      J = J + JUMP2 
 1010 CONTINUE
*
      IF (MOD(N,2).EQ.0) THEN 
      I = (N/2)*INC1 + 1
      J = (N/2)*INC2 + 1
*
      IF (IFORM.EQ.+1) THEN
      J = N*INC2 + 1
      JJ = J + INC2 
      DO 1015 L = 1 , LOT
      B(JJ) = 0.0
      JJ = JJ + JUMP2
 1015 CONTINUE
      ENDIF
*
      DO 1020 L = 1 , LOT
      B(J) = S * A(I)
      I = I + JUMP1 
      J = J + JUMP2 
 1020 CONTINUE
      ENDIF
*
*     FOLD ALONG LINES
*     ----------------
      ISHIFT = 2**24
      IQ = ILIST(1) / ISHIFT
      II = 3
      DO 1150 K = 1 , IQ
      IA = ILIST(II)
      JA = ILIST(II+1)
      IB = ILIST(II+2)
      JB = ILIST(II+3)
      S1 = SIGN(S,FLOAT(IB-IA))
      I = 1
      J = 1
      LV = LV1
      DO 1120 NB = 1 , NBLOX
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1100 L = 1 , LV
      V1(L) = S * A(IA+I)
      V2(L) = S1 * A(IB+I)
      I = I + JUMP1 
 1100 CONTINUE
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1110 L = 1 , LV
      B(JA+J) = V1(L)
      B(JB+J) = V2(L)
      J = J + JUMP2 
 1110 CONTINUE
      LV = MAXLOT
 1120 CONTINUE
      II = II + 4
 1150 CONTINUE
*
*     FOLD ON PLANES
*     --------------
      IQ = MOD(ILIST(1),ISHIFT)
      IF (IQ.EQ.0) GO TO 1500 
      DO 1250 K = 1 , IQ
      IA = ILIST(II)
      JA = ILIST(II+1)
      IB = ILIST(II+2)
      JB = ILIST(II+3)
      IC = ILIST(II+4)
      JC = ILIST(II+5)
      ID = ILIST(II+6)
      JD = ILIST(II+7)
      S1 = SIGN(S,FLOAT(ID-IA))
      S2 = SIGN(S,FLOAT(IC-IB))
      I = 1
      J = 1
      LV = LV1
      DO 1220 NB = 1 , NBLOX
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1200 L = 1 , LV
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      V1(L) = S  * ( AA + AD ) 
      V2(L) = S  * ( AA - AD ) 
      V3(L) = S2 * ( AC - AB )
      V4(L) = S1 * ( AC + AB )
      I = I + JUMP1 
 1200 CONTINUE
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1210 L = 1 , LV
      B(JA+J) = V2(L)
      B(JB+J) = V1(L)
      B(JC+J) = V3(L)
      B(JD+J) = V4(L)
      J = J + JUMP2 
 1210 CONTINUE
      LV = MAXLOT
 1220 CONTINUE
      II = II + 8
 1250 CONTINUE
*
*     FOLD ON CUBOIDS
*     ---------------
      IQ = ILIST(2)/ISHIFT
      IF (IQ.EQ.0) GO TO 1500 
      DO 1350 K = 1 , IQ
      IA = ILIST(II)
      JA = ILIST(II+1)
      IB = ILIST(II+2)
      JB = ILIST(II+3)
      IC = ILIST(II+4)
      JC = ILIST(II+5)
      ID = ILIST(II+6)
      JD = ILIST(II+7)
      IE = ILIST(II+8)
      JE = ILIST(II+9)
      IF = ILIST(II+10)
      JF = ILIST(II+11)
      IG = ILIST(II+12)
      JG = ILIST(II+13)
      IH = ILIST(II+14)
      JH = ILIST(II+15)
      S1 = SIGN(S,FLOAT(IH-IA))
      S2 = SIGN(S,FLOAT(IG-IB))
      S3 = SIGN(S,FLOAT(IF-IC))
      S4 = SIGN(S,FLOAT(IE-ID))
      I = 1
      J = 1
      LV = LV1
      DO 1320 NB = 1 , NBLOX
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1300 L = 1 , LV
      AA = A(IA+I)
      AB = A(IB+I)
      AC = A(IC+I)
      AD = A(ID+I)
      AE = A(IE+I)
      AF = A(IF+I)
      AG = A(IG+I)
      AH = A(IH+I)
      V1(L) = AA + AD
      V2(L) = AA - AD
      V3(L) = AC + AB
      V4(L) = AC - AB
      V5(L) = AE + AH
      V6(L) = AE - AH
      V7(L) = AG + AF
      V8(L) = AG - AF
      I = I + JUMP1 
 1300 CONTINUE
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1310 L = 1 , LV
      B(JA+J) = S * ( V2(L) - V7(L) )
      B(JD+J) = S * ( V2(L) + V7(L) )
      B(JB+J) = S * ( V1(L) - V8(L) )
      B(JC+J) = S * ( V1(L) + V8(L) )
      B(JE+J) = S4 * ( V6(L) - V3(L) )
      B(JH+J) = S1 * ( V6(L) + V3(L) )
      B(JF+J) = S3 * ( V5(L) - V4(L) )
      B(JG+J) = S2 * ( V5(L) + V4(L) )
      J = J + JUMP2 
 1310 CONTINUE
      LV = MAXLOT
 1320 CONTINUE
      II = II + 16
 1350 CONTINUE
*
*     FOLD ON HYPERCUBOIDS
*     --------------------
      IQ = MOD(ILIST(2),ISHIFT)
      IX = II
      IF (IQ.EQ.0) GO TO 1500 
      DO 1450 K = 1 , IQ
      IA = ILIST(IX)
      JA = ILIST(IX+1)
      IB = ILIST(IX+2)
      JB = ILIST(IX+3)
      IC = ILIST(IX+4)
      JC = ILIST(IX+5)
      ID = ILIST(IX+6)
      JD = ILIST(IX+7)
      IE = ILIST(IX+8)
      JE = ILIST(IX+9)
      IF = ILIST(IX+10)
      JF = ILIST(IX+11)
      IG = ILIST(IX+12)
      JG = ILIST(IX+13)
      IH = ILIST(IX+14)
      JH = ILIST(IX+15)
      II = ILIST(IX+16)
      JI = ILIST(IX+17)
      IJ = ILIST(IX+18)
      JJ = ILIST(IX+19)
      IK = ILIST(IX+20)
      JK = ILIST(IX+21)
      IL = ILIST(IX+22)
      JL = ILIST(IX+23)
      IM = ILIST(IX+24)
      JM = ILIST(IX+25)
      IN = ILIST(IX+26)
      JN = ILIST(IX+27)
      IO = ILIST(IX+28)
      JO = ILIST(IX+29)
      IP = ILIST(IX+30)
      JP = ILIST(IX+31)
      S1 = SIGN(S,FLOAT(IP-IA))
      S2 = SIGN(S,FLOAT(IO-IB))
      S3 = SIGN(S,FLOAT(IN-IC))
      S4 = SIGN(S,FLOAT(IM-ID))
      S5 = SIGN(S,FLOAT(IL-IE))
      S6 = SIGN(S,FLOAT(IK-IF))
      S7 = SIGN(S,FLOAT(IJ-IG))
      S8 = SIGN(S,FLOAT(II-IH))
*
      I = 1
      J = 1
      LV = LV1
      DO 1420 NB = 1 , NBLOX
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1400 L = 1 , LV
      AA = A(IA+I)
      AD = A(ID+I)
      AF = A(IF+I)
      AG = A(IG+I)
      T1 = AA + AD
      T2 = AA - AD
      T3 = AG + AF
      T4 = AG - AF
      V1(L) = S * ( T2 - T3 ) 
      V2(L) = S * ( T2 + T3 ) 
      V3(L) = S * ( T1 - T4 ) 
      V4(L) = S * ( T1 + T4 ) 
      AJ = A(IJ+I)
      AK = A(IK+I)
      AM = A(IM+I)
      AP = A(IP+I)
      T1 = AK + AJ
      T2 = AK - AJ
      T3 = AM + AP
      T4 = AM - AP
      V5(L) = S * ( T1 + T4 ) 
      V6(L) = S * ( T4 - T1 ) 
      V7(L) = S * ( T2 + T3 ) 
      V8(L) = S * ( T3 - T2 ) 
      AB = A(IB+I)
      AC = A(IC+I)
      AE = A(IE+I)
      AH = A(IH+I)
      AI = A(II+I)
      AL = A(IL+I)
      AN = A(IN+I)
      AO = A(IO+I)
      V9(L)  = AC + AB
      V10(L) = AC - AB
      V11(L) = AE + AH
      V12(L) = AE - AH
      V13(L) = AI + AL
      V14(L) = AI - AL
      V15(L) = AO + AN
      V16(L) = AO - AN
      I = I + JUMP1 
 1400 CONTINUE
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1410 L = 1 , LV
      B(JA+J) = V1(L) - V5(L) 
      B(JH+J) = V1(L) + V5(L) 
      B(JC+J) = V4(L) - V8(L) 
      B(JF+J) = V4(L) + V8(L) 
      B(JB+J) = V3(L) - V7(L) 
      B(JG+J) = V3(L) + V7(L) 
      B(JD+J) = V2(L) - V6(L) 
      B(JE+J) = V2(L) + V6(L) 
      T1 = V12(L) + V9(L)
      T2 = V12(L) - V9(L)
      T3 = V14(L) - V15(L)
      T4 = V14(L) + V15(L)
      B(JI+J) = S8 * ( T3 - T1 )
      B(JP+J) = S1 * ( T3 + T1 )
      B(JL+J) = S5 * ( T4 - T2 )
      B(JM+J) = S4 * ( T4 + T2 )
      T1 = V11(L) + V10(L)
      T2 = V11(L) - V10(L)
      T3 = V13(L) - V16(L)
      T4 = V13(L) + V16(L)
      B(JJ+J) = S7 * ( T3 - T1 )
      B(JO+J) = S2 * ( T3 + T1 )
      B(JK+J) = S6 * ( T4 - T2 )
      B(JN+J) = S3 * ( T4 + T2 )
      J = J + JUMP2 
 1410 CONTINUE
      LV = MAXLOT
 1420 CONTINUE
      IX = IX + 32
 1450 CONTINUE
*
 1500 CONTINUE
      IERR = 0
 2000 CONTINUE
      RETURN
      END 
*DECK FORQFA
*        SUBROUTINE 'QPFA2'
*        SELF-SORTING IN-PLACE PRIME FACTOR FFT
*         (HALF-COMPLEX TO REAL)
*
*        CALL QPFA2(A,B,IFAX,ILIST,INC1,JUMP1,INC2,JUMP2,
*       +           N,LOT,IFORM,IERR)
*
*        A IS THE INPUT ARRAY (SPECTRAL COEFFICIENTS)
*        B IS THE OUTPUT ARRAY (GRIDPOINT VALUES) 
*        IFAX IS A LIST OF FACTORS OF N 
*              IFAX(1) CONTAINS THE NUMBER OF FACTORS (NFAX)
*              THE FACTORS FOLLOW IN IFAX(2) TO IFAX(NFAX+1)
*        ILIST IS A LIST OF ADDRESSES USED IN THE 'FOLDING' 
*          STAGE OF THE ALGORITHM
*        INC1 IS THE INCREMENT WITHIN THE INPUT DATA VECTORS
*        JUMP1 IS THE INCREMENT BETWEEN THE INPUT DATA VECTORS
*        INC2 IS THE INCREMENT WITHIN THE OUTPUT DATA VECTORS
*        JUMP2 IS THE INCREMENT BETWEEN THE OUTPUT DATA VECTORS
*        N IS THE LENGTH OF THE TRANSFORMS
*        LOT IS THE NUMBER OF TRANSFORMS
*        IFORM = +1 FOR CONVENTIONAL ORDERING OF SPECTRAL COEFFS
*                -1 FOR UNCONVENTIONAL (FAST-IN-PLACE) ORDERING
*        IERR IS AN ERROR INDICATOR: 
*                  0 - TRANSFORM COMPLETED WITHOUT ERROR
*                  1 - FACTORS ARE NOT MUTUALLY PRIME
*                  2 - ILLEGAL FACTOR
*
*        FORTRAN VERSION WRITTEN BY CLIVE TEMPERTON
*        RECHERCHE EN PREVISION NUMERIQUE -- NOVEMBER 1986
*
*        MODIFIED BY CLIVE TEMPERTON MAY 1987 -- TO GUARANTEE
*        CORRECT IN-PLACE WORKING OF FORTRAN VERSION
*
*        MODIFIED BY BERNARD DUGAS DECEMBER 1992 -- USE "IMPLICIT NONE"
*        AND USE VECTOR TEMPORARIES WITHIN LOOPS
*
*----------------------------------------------------------------------
*
*        FOR A MATHEMATICAL DEVELOPMENT OF THE ALGORITHM USED,
*        SEE: 
*
*        C TEMPERTON : "IMPLEMENTATION OF A SELF-SORTING IN-PLACE
*          PRIME FACTOR FFT ALGORITHM", 
*          JOURNAL OF COMPUTATIONAL PHYSICS VOL 58 (1985), 283-299.
*
*        C TEMPERTON : "A SELF-SORTING IN-PLACE PRIME FACTOR
*          REAL/HALF-COMPLEX FFT ALGORITHM",
*          TO APPEAR IN JOURNAL OF COMPUTATIONAL PHYSICS
*
*        C TEMPERTON : "A NEW SET OF MINIMUM-ADD SMALL-N DFT MODULES",
*          TO APPEAR IN JOURNAL OF COMPUTATIONAL PHYSICS
*
*----------------------------------------------------------------------
*
      SUBROUTINE QPFA2(A,B,IFAX,ILIST,INC1,JUMP1,INC2,JUMP2, 
     *                 N,LOT,IFORM,IERR) 
*
      IMPLICIT    none
*
      INTEGER     INC1,JUMP1,INC2,JUMP2,N,LOT,IFORM,IERR
      REAL        A(*),B(*)
      INTEGER     IFAX(*),ILIST(*)
*
      INTEGER     MAXLOT
      PARAMETER ( MAXLOT = lot_maximum )

      REAL       V1(MAXLOT+1), V2(MAXLOT+1), V3(MAXLOT+1), V4(MAXLOT+1),
     *           V5(MAXLOT+1), V6(MAXLOT+1), V7(MAXLOT+1), V8(MAXLOT+1),
     *           V9(MAXLOT+1),V10(MAXLOT+1),V11(MAXLOT+1),V12(MAXLOT+1),
     *          V13(MAXLOT+1),V14(MAXLOT+1),V15(MAXLOT+1),V16(MAXLOT+1)
*
      REAL        AA,AB,AC,AD,AE,AF,AG,AH,AI,AJ,AK,AL,AM,AN,AO,AP
      REAL        BA,BB,BC,BD,BE,BF,BG,BH,BI,BJ,BK,BL,   BN,BO,BP
      REAL        S,S1,S2,S3,S4,S5,S6,S7,S8,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9
      REAL        T1,T2,T3,T4,T5,T6
*
      INTEGER     NBLOX,I,J,K,L,LV,LV1,M,MM,MU,NN,NB
      INTEGER     IA,IB,IC,ID,IE,IF,IG,IH,II,IJ,IK,IL,IM,IN,IO,IP,IQ
      INTEGER     JA,JB,JC,JD,JE,JF,JG,JH,JI,JJ,JK,JL,JM,JN,JO,JP
      INTEGER     NFAX,IFAC,ISHIFT,IX
*
      REAL   SIN60,QRT5, SIN36,SIN72,SIN45
      REAL   COS26,SIN26,COS51,SIN51
      REAL   COS77,SIN77,COS20,SIN20
      REAL   COS40,SIN40,COS80,SIN80
      REAL   COS23,SIN23
*
      DATA   SIN60/0.866025403784437/
      DATA   QRT5 /0.559016994374947/
      DATA   SIN36/0.587785252292473/,SIN72/0.951056516295154/
      DATA   COS26/0.900968867902419/,SIN26/0.433883739117558/
      DATA   COS51/0.623489801858734/,SIN51/0.781831482468030/
      DATA   COS77/0.222520933956314/,SIN77/0.974927912181824/
      DATA   SIN45/0.707106781186548/
      DATA   SIN20/0.342020143325669/,COS20/0.939692620785908/,
     *       SIN40/0.642787609686539/,COS40/0.766044443118978/,
     *       SIN80/0.984807753012208/,COS80/0.173648177666930/
      DATA   SIN23/0.382683432365090/,COS23/0.923879532511287/
*
      NBLOX = 1 + (LOT-1)/MAXLOT
      LV1 = LOT - (NBLOX-1)*MAXLOT
*
*     K=0 AND K=N/2 
*     ------------- 
      I = 1
      J = 1
      DO 1010 L = 1 , LOT
      B(I) = A(J)
      I = I + JUMP2 
      J = J + JUMP1 
 1010 CONTINUE
      IF (MOD(N,2).EQ.0) THEN 
      I = (N/2)*INC2 + 1
      J = (N/2)*INC1 + 1
      IF (IFORM.EQ.+1) J = N*INC1 + 1
      DO 1020 L = 1 , LOT
      B(I) = A(J)
      I = I + JUMP2 
      J = J + JUMP1 
 1020 CONTINUE
      ENDIF
*
*     FOLD ALONG LINES
*     ----------------
      S = 2.0
      ISHIFT = 2**24
      IQ = ILIST(1)/ISHIFT
      II = 3
      DO 1150 K = 1 , IQ
      IA = ILIST(II)
      JA = ILIST(II+1)
      IB = ILIST(II+2)
      JB = ILIST(II+3)
      S1 = SIGN(S,FLOAT(IB-IA))
      I = 1
      J = 1
      LV = LV1
      DO 1120 NB = 1 , NBLOX
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1100 L = 1 , LV
      V1(L) = S * A(JA+J)
      V2(L) = S1 * A(JB+J)
      J = J + JUMP1 
 1100 CONTINUE
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1110 L = 1 , LV
      B(IA+I) = V1(L)
      B(IB+I) = V2(L)
      I = I + JUMP2 
 1110 CONTINUE
      LV = MAXLOT
 1120 CONTINUE
      II = II + 4
 1150 CONTINUE
*
*     FOLD ON PLANES
*     --------------
      IQ = MOD(ILIST(1),ISHIFT)
      IF (IQ.EQ.0) GO TO 1500 
      DO 1250 K = 1 , IQ
      IA = ILIST(II)
      JA = ILIST(II+1)
      IB = ILIST(II+2)
      JB = ILIST(II+3)
      IC = ILIST(II+4)
      JC = ILIST(II+5)
      ID = ILIST(II+6)
      JD = ILIST(II+7)
      S1 = SIGN(S,FLOAT(ID-IA))
      S2 = SIGN(S,FLOAT(IC-IB))
      I = 1
      J = 1
      LV = LV1
      DO 1220 NB = 1 , NBLOX
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1200 L = 1 , LV
      AA = A(JA+J)
      AB = A(JB+J)
      AC = A(JC+J)
      AD = A(JD+J)
      V1(L) = S  * ( AB + AA ) 
      V2(L) = S  * ( AB - AA ) 
      V3(L) = S1 * AD - S2 * AC
      V4(L) = S1 * AD + S2 * AC
      J = J + JUMP1 
 1200 CONTINUE
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1210 L = 1 , LV
      B(IA+I) = V1(L)
      B(ID+I) = V2(L)
      B(IB+I) = V3(L)
      B(IC+I) = V4(L)
      I = I + JUMP2 
 1210 CONTINUE
      LV = MAXLOT
 1220 CONTINUE
      II = II + 8
 1250 CONTINUE
*
*     FOLD ON CUBOIDS
*     ---------------
      IQ = ILIST(2)/ISHIFT
      IF (IQ.EQ.0) GO TO 1500 
      DO 1350 K = 1 , IQ
      IA = ILIST(II)
      JA = ILIST(II+1)
      IB = ILIST(II+2)
      JB = ILIST(II+3)
      IC = ILIST(II+4)
      JC = ILIST(II+5)
      ID = ILIST(II+6)
      JD = ILIST(II+7)
      IE = ILIST(II+8)
      JE = ILIST(II+9)
      IF = ILIST(II+10)
      JF = ILIST(II+11)
      IG = ILIST(II+12)
      JG = ILIST(II+13)
      IH = ILIST(II+14)
      JH = ILIST(II+15)
      S1 = SIGN(S,FLOAT(IH-IA))
      S2 = SIGN(S,FLOAT(IG-IB))
      S3 = SIGN(S,FLOAT(IF-IC))
      S4 = SIGN(S,FLOAT(IE-ID))
      I = 1
      J = 1
      LV = LV1
      DO 1320 NB = 1 , NBLOX
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1300 L = 1 , LV
      AA = A(JA+J)
      AB = A(JB+J)
      AC = A(JC+J)
      AD = A(JD+J)
      AE = A(JE+J)
      AF = A(JF+J)
      AG = A(JG+J)
      AH = A(JH+J)
      V1(L) = S  * ( AD + AA ) 
      V8(L) = S  * ( AD - AA ) 
      V4(L) = S1 * AH - S4 * AE
      V5(L) = S1 * AH + S4 * AE 
      V2(L) = S  * ( AC + AB ) 
      V7(L) = S  * ( AC - AB ) 
      V3(L) = S2 * AG - S3 * AF
      V6(L) = S2 * AG + S3 * AF
      J = J + JUMP1 
 1300 CONTINUE
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1310 L = 1 , LV
      B(IA+I) = V2(L) + V1(L) 
      B(ID+I) = V2(L) - V1(L) 
      B(IB+I) = V4(L) - V3(L) 
      B(IC+I) = V4(L) + V3(L) 
      B(IE+I) = V6(L) + V5(L) 
      B(IH+I) = V6(L) - V5(L) 
      B(IF+I) = V8(L) - V7(L) 
      B(IG+I) = V8(L) + V7(L) 
      I = I + JUMP2 
 1310 CONTINUE
      LV = MAXLOT
 1320 CONTINUE
      II = II + 16
 1350 CONTINUE
*
*     FOLD ON HYPERCUBOIDS
*     --------------------
      IQ = MOD(ILIST(2),ISHIFT)
      IX = II
      IF (IQ.EQ.0) GO TO 1500 
      DO 1450 K = 1 , IQ
      IA = ILIST(IX)
      JA = ILIST(IX+1)
      IB = ILIST(IX+2)
      JB = ILIST(IX+3)
      IC = ILIST(IX+4)
      JC = ILIST(IX+5)
      ID = ILIST(IX+6)
      JD = ILIST(IX+7)
      IE = ILIST(IX+8)
      JE = ILIST(IX+9)
      IF = ILIST(IX+10)
      JF = ILIST(IX+11)
      IG = ILIST(IX+12)
      JG = ILIST(IX+13)
      IH = ILIST(IX+14)
      JH = ILIST(IX+15)
      II = ILIST(IX+16)
      JI = ILIST(IX+17)
      IJ = ILIST(IX+18)
      JJ = ILIST(IX+19)
      IK = ILIST(IX+20)
      JK = ILIST(IX+21)
      IL = ILIST(IX+22)
      JL = ILIST(IX+23)
      IM = ILIST(IX+24)
      JM = ILIST(IX+25)
      IN = ILIST(IX+26)
      JN = ILIST(IX+27)
      IO = ILIST(IX+28)
      JO = ILIST(IX+29)
      IP = ILIST(IX+30)
      JP = ILIST(IX+31)
      S1 = SIGN(S,FLOAT(IP-IA))
      S2 = SIGN(S,FLOAT(IO-IB))
      S3 = SIGN(S,FLOAT(IN-IC))
      S4 = SIGN(S,FLOAT(IM-ID))
      S5 = SIGN(S,FLOAT(IL-IE))
      S6 = SIGN(S,FLOAT(IK-IF))
      S7 = SIGN(S,FLOAT(IJ-IG))
      S8 = SIGN(S,FLOAT(II-IH))
*
      I = 1
      J = 1
      LV = LV1
      DO 1420 NB = 1 , NBLOX
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1400 L = 1 , LV
      AA = A(JA+J)
      AD = A(JD+J)
      AE = A(JE+J)
      AH = A(JH+J)
      T1 = AH + AA
      T2 = AH - AA
      T3 = AE + AD
      T4 = AE - AD
      V1(L) = S * ( T1 + T3 ) 
      V8(L) = S * ( T3 - T1 ) 
      V4(L) = S * ( T2 - T4 ) 
      V5(L) = S * ( T2 + T4 ) 
      AI = A(JI+J)
      AL = A(JL+J)
      AM = A(JM+J)
      AP = A(JP+J)
      T1 = S1*AP - S8*AI
      T2 = S1*AP + S8*AI
      T3 = S4*AM - S5*AL
      T4 = S4*AM + S5*AL
      V9(L)  = T1 - T3
      V12(L) = T1 + T3
      V13(L) = T4 + T2
      V16(L) = T4 - T2
      AB = A(JB+J)
      AC = A(JC+J)
      AF = A(JF+J)
      AG = A(JG+J)
      T1 = AG + AB
      T2 = AG - AB
      T3 = AF + AC
      T4 = AF - AC
      V2(L) = S * ( T1 + T3 ) 
      V7(L) = S * ( T3 - T1 ) 
      V3(L) = S * ( T2 - T4 ) 
      V6(L) = S * ( T2 + T4 ) 
      AJ = A(JJ+J)
      AK = A(JK+J)
      AN = A(JN+J)
      AO = A(JO+J)
      T1 = S2*AO - S7*AJ
      T2 = S2*AO + S7*AJ
      T3 = S3*AN - S6*AK
      T4 = S3*AN + S6*AK
      V10(L) = T1 - T3
      V11(L) = T1 + T3
      V14(L) = T4 + T2
      V15(L) = T4 - T2
      J = J + JUMP1 
 1400 CONTINUE
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1410 L = 1 , LV
      B(IA+I) = V2(L) + V1(L) 
      B(ID+I) = V2(L) - V1(L) 
      B(IJ+I) = V4(L) - V3(L) 
      B(IK+I) = V4(L) + V3(L) 
      B(IB+I) = V9(L) - V10(L)
      B(IC+I) = V9(L) + V10(L)
      B(II+I) = V14(L) + V13(L)
      B(IL+I) = V14(L) - V13(L)
      B(IM+I) = V6(L) + V5(L) 
      B(IP+I) = V6(L) - V5(L) 
      B(IE+I) = V11(L) + V12(L)
      B(IH+I) = V11(L) - V12(L)
      B(IF+I) = V8(L) - V7(L) 
      B(IG+I) = V8(L) + V7(L) 
      B(IN+I) = V16(L) - V15(L)
      B(IO+I) = V16(L) + V15(L)
      I = I + JUMP2 
 1410 CONTINUE
      LV = MAXLOT
 1420 CONTINUE
      IX = IX + 32
 1450 CONTINUE
*
 1500 CONTINUE
C
      NFAX = IFAX(1)
*
      DO 1000 K = 1 , NFAX
C
      IFAC=IFAX(K+1)
      M=N/IFAC
      DO 100 J = 1 , IFAC
      MU=J
      MM=J*M
      IF (MOD(MM,IFAC).EQ.1) GO TO 110
  100 CONTINUE
*
*     ERROR EXIT (FACTORS NOT MUTUALLY PRIME)
*     ----------
      IERR = 1
      GO TO 2000
*
  110 CONTINUE
      MM = MM * INC2
      NN = N * INC2 
C
      IA=1
      IB=IA+MM
      IF (IFAC.EQ.2) GO TO 200
      IC=IB+MM
      IF (IC.GT.NN) IC=IC-NN
      IF (IFAC.EQ.3) GO TO 300
      ID=IC+MM
      IF (ID.GT.NN) ID=ID-NN
      IF (IFAC.EQ.4) GO TO 400
      IE=ID+MM
      IF (IE.GT.NN) IE=IE-NN
      IF (IFAC.EQ.5) GO TO 500
      IF=IE+MM
      IF (IF.GT.NN) IF=IF-NN
      IG=IF+MM
      IF (IG.GT.NN) IG=IG-NN
      IF (IFAC.EQ.7) GO TO 700
      IH=IG+MM
      IF (IH.GT.NN) IH=IH-NN
      IF (IFAC.EQ.8) GO TO 800
      II=IH+MM
      IF (II.GT.NN) II=II-NN
      IF (IFAC.EQ.9) GO TO 900
      IJ=II+MM
      IF (IJ.GT.NN) IJ=IJ-NN
      IK=IJ+MM
      IF (IK.GT.NN) IK=IK-NN
      IL=IK+MM
      IF (IL.GT.NN) IL=IL-NN
      IM=IL+MM
      IF (IM.GT.NN) IM=IM-NN
      IN=IM+MM
      IF (IN.GT.NN) IN=IN-NN
      IO=IN+MM
      IF (IO.GT.NN) IO=IO-NN
      IP=IO+MM
      IF (IP.GT.NN) IP=IP-NN
      IF (IFAC.EQ.16) GO TO 1600
*
*     ERROR EXIT (ILLEGAL FACTOR)
*     ----------
      IERR = 2
      GO TO 2000
C
C     CODING FOR FACTOR 2
C     -------------------
  200 CONTINUE
C
      DO 220 L = 1 , M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 210 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      B(IB+I) = BA - BB
      B(IA+I) = BA + BB
      I = I + JUMP2 
  210 CONTINUE
      IX=IB+INC2
      IB=IA+INC2
      IA=IX
  220 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 3
C     -------------------
  300 CONTINUE
      Z3=SIN60
      IF (MU.EQ.2) Z3=-Z3
C
      DO 320 L = 1 , M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 310 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      T1 = BA - .5 * BB
      T2 = Z3 * BC
      B(IA+I) = BA + BB
      B(IB+I) = T1 - T2
      B(IC+I) = T1 + T2
      I = I + JUMP2 
  310 CONTINUE
      IX=IC+INC2
      IC=IB+INC2
      IB=IA+INC2
      IA=IX
  320 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 4
C     -------------------
  400 CONTINUE
      Z4=1.0
      IF (MU.EQ.3) Z4=-Z4
C
      DO 420 L = 1 , M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 410 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      BD = B(ID+I)
      T1 = BA + BC
      T2 = BA - BC
      B(IA+I) = T1 + BB
      B(IC+I) = T1 - BB
      B(IB+I) = T2 - Z4 * BD
      B(ID+I) = T2 + Z4 * BD
      I = I + JUMP2 
  410 CONTINUE
      IX=ID+INC2
      ID=IC+INC2
      IC=IB+INC2
      IB=IA+INC2
      IA=IX
  420 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 5
C     -------------------
  500 CONTINUE
      IF (MU.EQ.1) THEN
        Z1=QRT5
        Z2=SIN72
        Z3=SIN36
      ELSE IF (MU.EQ.2) THEN
        Z1=-QRT5
        Z2=SIN36
        Z3=-SIN72
      ELSE IF (MU.EQ.3) THEN
        Z1=-QRT5
        Z2=-SIN36
        Z3=SIN72
      ELSE IF (MU.EQ.4) THEN
        Z1=QRT5
        Z2=-SIN72
        Z3=-SIN36
      ENDIF
      DO 520 L=1,M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 510 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      BD = B(ID+I)
      BE = B(IE+I)
      T1 = BB + BC
      T2 = Z1 * ( BB - BC )
      T3 = BA - 0.25 * T1
      B(IA+I) =   BA + T1
      T4 = T3 + T2
      T5 = T3 - T2
      T1 = Z2 * BE + Z3 * BD
      T2 = Z3 * BE - Z2 * BD
      B(IB+I) = T4 - T1
      B(IE+I) = T4 + T1
      B(IC+I) = T5 - T2
      B(ID+I) = T5 + T2
      I = I + JUMP2 
  510 CONTINUE
      IX = IE + INC2
      IE = ID + INC2
      ID = IC + INC2
      IC = IB + INC2
      IB = IA + INC2
      IA = IX
  520 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 7
C     -------------------
  700 CONTINUE
      IF (MU.EQ.1) THEN
        Z1=COS51
        Z2=-COS77
        Z3=-COS26
        Z4=SIN51
        Z5=SIN77
        Z6=SIN26
      ELSE IF (MU.EQ.2) THEN
        Z1=-COS77
        Z2=-COS26
        Z3=COS51
        Z4=SIN77
        Z5=-SIN26
        Z6=-SIN51
      ELSE IF (MU.EQ.3) THEN
        Z1=-COS26
        Z2=COS51
        Z3=-COS77
        Z4=SIN26
        Z5=-SIN51
        Z6=SIN77
      ELSE IF (MU.EQ.4) THEN
        Z1=-COS26
        Z2=COS51
        Z3=-COS77
        Z4=-SIN26
        Z5=SIN51
        Z6=-SIN77
      ELSE IF (MU.EQ.5) THEN
        Z1=-COS77
        Z2=-COS26
        Z3=COS51
        Z4=-SIN77
        Z5=SIN26
        Z6=SIN51
      ELSE IF (MU.EQ.6) THEN
        Z1=COS51
        Z2=-COS77
        Z3=-COS26
        Z4=-SIN51
        Z5=-SIN77
        Z6=-SIN26
      ENDIF
      DO 720 L=1,M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 710 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      BD = B(ID+I)
      BE = B(IE+I)
      BF = B(IF+I)
      BG = B(IG+I)
      T1 = BA + Z1 * BB + Z2 * BC + Z3 * BD
      T2 = BA + Z2 * BB + Z3 * BC + Z1 * BD
      T3 = BA + Z3 * BB + Z1 * BC + Z2 * BD
      T4 =      Z4 * BG + Z5 * BF + Z6 * BE
      T5 =      Z5 * BG - Z6 * BF - Z4 * BE
      T6 =      Z6 * BG - Z4 * BF + Z5 * BE
      B(IA+I) = BA + BB + BC + BD
      B(IB+I) = T1 - T4
      B(IC+I) = T2 - T5
      B(ID+I) = T3 - T6
      B(IE+I) = T3 + T6
      B(IF+I) = T2 + T5
      B(IG+I) = T1 + T4
      I = I + JUMP2 
  710 CONTINUE
      IX = IG + INC2
      IG = IF + INC2
      IF = IE + INC2
      IE = ID + INC2
      ID = IC + INC2
      IC = IB + INC2
      IB = IA + INC2
      IA = IX
  720 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 8
C     -------------------
  800 CONTINUE
      IF (MU.EQ.1) THEN
        Z1=1.0
        Z2=SIN45
        Z3=SIN45
      ELSE IF (MU.EQ.3) THEN
        Z1=-1.0
        Z2=-SIN45
        Z3=SIN45
      ELSE IF (MU.EQ.5) THEN
        Z1=1.0
        Z2=-SIN45
        Z3=-SIN45
      ELSE IF (MU.EQ.7) THEN
        Z1=-1.0
        Z2=SIN45
        Z3=-SIN45
      ENDIF
      DO 820 L=1,M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 810 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      BD = B(ID+I)
      BE = B(IE+I)
      BF = B(IF+I)
      BG = B(IG+I)
      BH = B(IH+I)
      T1 = BA + BE
      T2 = BA - BE
      T3 = T1 + BC
      T4 = T1 - BC
      T5 = BB + BD
      T6 = Z2 * ( BB - BD )
      B(IA+I) = T3 + T5
      B(IE+I) = T3 - T5
      T1 = T2 + T6
      T3 = T2 - T6
      T5 = Z1 * ( BH - BF )
      T6 = Z3 * ( BH + BF )
      T2 = Z1 * B(IG+I)
      B(IC+I) = T4 - T5
      B(IG+I) = T4 + T5
      T4 = T6 - T2
      T5 = T6 + T2
      B(ID+I) = T3 - T4
      B(IF+I) = T3 + T4
      B(IB+I) = T1 - T5
      B(IH+I) = T1 + T5
      I = I + JUMP2 
  810 CONTINUE
      IX = IH + INC2
      IH = IG + INC2
      IG = IF + INC2
      IF = IE + INC2
      IE = ID + INC2
      ID = IC + INC2
      IC = IB + INC2
      IB = IA + INC2
      IA = IX
  820 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 9
C     -------------------
  900 CONTINUE
      IF (MU.EQ.1) THEN
        Z1=SIN60
        Z2=COS40
        Z3=SIN40
        Z4=COS80
        Z5=SIN80
      ELSE IF (MU.EQ.2) THEN
        Z1=-SIN60
        Z2=COS80
        Z3=SIN80
        Z4=-COS20
        Z5=SIN20
      ELSE IF (MU.EQ.4) THEN
        Z1=SIN60
        Z2=-COS20
        Z3=SIN20
        Z4=COS40
        Z5=-SIN40
      ELSE IF (MU.EQ.5) THEN
        Z1=-SIN60
        Z2=-COS20
        Z3=-SIN20
        Z4=COS40
        Z5=SIN40
      ELSE IF (MU.EQ.7) THEN
        Z1=SIN60
        Z2=COS80
        Z3=-SIN80
        Z4=-COS20
        Z5=-SIN20
      ELSE IF (MU.EQ.8) THEN
        Z1=-SIN60
        Z2=COS40
        Z3=-SIN40
        Z4=COS80
        Z5=-SIN80
      ENDIF
      Z6=Z1*Z2
      Z7=Z1*Z3
      Z8=Z1*Z4
      Z9=Z1*Z5
      DO 920 L = 1 , M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 910 J = 1 , LOT
      BA = B(IA+I)
      BB = B(IB+I)
      BC = B(IC+I)
      BD = B(ID+I)
      BE = B(IE+I)
      BF = B(IF+I)
      BG = B(IG+I)
      BH = B(IH+I)
      BI = B(II+I)
      T1 = BA - .5 * BD
      BA = BA + BD
      T4 = T1 - Z1 * BG
      T2 = T1 + Z1 * BG
      T1 = BE + BC
      BC = Z1 * ( BE - BC )
      BE = BB - .5 * T1
      BB = BB + T1
      T1 = BF - BH
      BF = Z1 * ( BF + BH )
      BH = BI - .5 * T1
      BI = Z1 * ( BI + T1 )
      T1 = BA - .5 * BB
      B(IA+I) = BA + BB
      B(ID+I) = T1 - BI
      B(IG+I) = T1 + BI
      BI = T4
      T3 = BE - BF
      T4 = BE + BF
      T5 = BH + BC
      T6 = BH - BC
      T1 = Z2 * T3 - Z3 * T5
      BF = Z7 * T3 + Z6 * T5
      B(IB+I) = BI + T1
      T1 = BI - .5 * T1 
      B(IE+I) = T1 - BF
      B(IH+I) = T1 + BF
      T3 = Z4 * T4 - Z5 * T6
      T5 = Z9 * T4 + Z8 * T6
      B(IC+I) = T2 + T3
      T1 = T2 - .5 * T3
      B(IF+I) = T1 - T5
      B(II+I) = T1 + T5
      I = I + JUMP2 
  910 CONTINUE
      IX = II + INC2
      II = IH + INC2
      IH = IG + INC2
      IG = IF + INC2
      IF = IE + INC2
      IE = ID + INC2
      ID = IC + INC2
      IC = IB + INC2
      IB = IA + INC2
      IA = IX
  920 CONTINUE
C
      GO TO 1000
C
C     CODING FOR FACTOR 16
C     --------------------
 1600 CONTINUE
      IF (MU.EQ.1) THEN
        Z1=1.0
        Z2=COS23
        Z3=SIN23
        Z4=SIN45
      ELSE IF (MU.EQ.3) THEN
        Z1=-1.0
        Z2=SIN23
        Z3=COS23
        Z4=-SIN45
      ELSE IF (MU.EQ.5) THEN
        Z1=1.0
        Z2=-SIN23
        Z3=COS23
        Z4=-SIN45
      ELSE IF (MU.EQ.7) THEN
        Z1=-1.0
        Z2=-COS23
        Z3=SIN23
        Z4=SIN45
      ELSE IF (MU.EQ.9) THEN
        Z1=1.0
        Z2=-COS23
        Z3=-SIN23
        Z4=SIN45
      ELSE IF (MU.EQ.11) THEN 
        Z1=-1.0
        Z2=-SIN23
        Z3=-COS23
        Z4=-SIN45
      ELSE IF (MU.EQ.13) THEN 
        Z1=1.0
        Z2=SIN23
        Z3=-COS23
        Z4=-SIN45
      ELSE IF (MU.EQ.15) THEN 
        Z1=-1.0
        Z2=COS23
        Z3=-SIN23
        Z4=SIN45
      ENDIF
      Z5=Z1*Z4
      Z6=Z1*Z3
      Z7=Z1*Z2
      DO 1620 L = 1 , M
      I = 0
#if defined (CRAY)
CDIR$ IVDEP
#endif
#if defined (NEC)
*VDIR NODEP
#endif
      DO 1610 J = 1 , LOT
      BA = B(IA+I)
      BE = B(IE+I)
      BI = B(II+I)
      T1 = BA + BI
      T2 = BA - BI
      BA = T1 + BE
      BE = T1 - BE
      BB = B(IB+I)
      BD = B(ID+I)
      BF = B(IF+I)
      BH = B(IH+I)
      T3 = BB + BH
      T4 = BB - BH
      T5 = BF + BD
      T6 = BF - BD
      BI = T3 + T5
      BF = Z4 * ( T3 - T5 )
      BB = Z2 * T4 - Z6 * T6
      BH = Z6 * T4 + Z2 * T6
      BC = B(IC+I)
      BG = B(IG+I)
      T3 = BC + BG
      T4 = Z4 * ( BC - BG )
      T5 = BA + T3
      T6 = BA - T3
      B(IA+I) = T5 + BI
      B(II+I) = T5 - BI
      T1 = T2 + T4
      T2 = T2 - T4
      BC = T1 - BB
      BB = T1 + BB
      BD = T2 + BH
      BH = T2 - BH
      BG = BE + BF
      BF = BE - BF
      BJ = B(IJ+I)
      BL = B(IL+I)
      BN = B(IN+I)
      BP = B(IP+I)
      T3 = BP - BJ
      T4 = BL - BN
      T5 = Z1 * ( T3 + T4 )
      T1 = Z5 * ( T3 - T4 )
      B(IE+I) = T6 - T5
      T2 = T6 + T5
      T3 = BP + BJ
      T4 = BL + BN
      BL = Z7 * T3 - Z3 * T4
      BN = Z3 * T3 + Z7 * T4
      BK = B(IK+I)
      BO = B(IO+I)
      T5 = Z1 * ( BO - BK )
      T6 = Z5 * ( BO + BK )
      T3 = T1 + T5
      T4 = T1 - T5
      T1 = Z1 * B(IM+I) - T6
      T5 = Z1 * B(IM+I) + T6
      B(IM+I) = T2
      T6 = BL + T1
      T1 = BL - T1
      T2 = BN - T5
      T5 = BN + T5
      B(IP+I) = BB + T5
      B(IB+I) = BB - T5
      B(IN+I) = BD + T1
      B(ID+I) = BD - T1
      B(IO+I) = BG + T3
      T1 = BG - T3
      B(IK+I) = BF + T4
      B(IG+I) = BF - T4
      B(IF+I) = BH - T6
      B(IL+I) = BH + T6
      B(IH+I) = BC - T2
      B(IJ+I) = BC + T2
      B(IC+I) = T1
      I = I + JUMP2 
 1610 CONTINUE
      IX = IP + INC2
      IP = IO + INC2
      IO = IN + INC2
      IN = IM + INC2
      IM = IL + INC2
      IL = IK + INC2
      IK = IJ + INC2
      IJ = II + INC2
      II = IH + INC2
      IH = IG + INC2
      IG = IF + INC2
      IF = IE + INC2
      IE = ID + INC2
      ID = IC + INC2
      IC = IB + INC2
      IB = IA + INC2
      IA = IX
 1620 CONTINUE
C
 1000 CONTINUE
      IERR = 0
 2000 CONTINUE
      RETURN
      END 

