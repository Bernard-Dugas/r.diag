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
C     $Log: ffwfg2.ftn,v $
C     Revision 3.2  2014/09/25 18:42:02  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.1  1994/11/17 14:13:19  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:55:29  13:55:29  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:31:42  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.0  92/02/21  11:32:46  armnrbd
C     Initial revision
C     

      SUBROUTINE ffwfg2 (FW,LFW,GR,LGR,IR,ILONG,WRKS,NL)

C     *****   FEB 1976  -  JOHN D. HENDERSON  ****

C     * MODIFIED NOV 24/82 - R.LAPRISE, S.J.LAMBERT.

C     * USE THE FORTRAN FAST FOURIER TRANSFORM TO PERFORM A GRID TO 
C     * WAVE FOURIER TRANSFORM OF MULTI-LEVEL FIELD.

C     * NOTE... GR AND FW ARE GENERALLY EQUIVALENCED IN CALLING ROUTINE

C     * PARAMETERS...
C     * FW    = MULTI-LEVEL COMPLEX ARRAY TO RETURN THE FOURIER 
C     *         AMPLITUDES
C     * LFW   = FIRST DIMENSION OF FW
C     * GR    = MULTI-LEVEL ARRAY CARRYING GRID POINT FIELDS.
C     * LGR   = FIRST DIMENSION OF GR.
C     * IR    = MAXIMUM E.W. WAVE NUMBER TO TRANSFORM TO.
C     * ILONG = NUMBER OF LONGITUDES ON THE TRANSFORM GRID.
C     * WRKS  = WORK FIELD.
C     * NL    = NUMBER OF VERTICAL LEVELS.

      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)

      COMPLEX    FW(LFW,NL)
      DIMENSION  GR(LGR,NL),WRKS(1)
      COMPLEX    AA,BB 

C-----------------------------------------------------------------------
      IRP1=IR+1
      IF (MOD(ILONG,3).NE.0) THEN
C-----------------------------------------------------------------------
C     * CASE 1 - ILONG IS A POWER OF 2. 

          IF (IR.GT.ILONG/2) THEN
              WRITE(6,6010) IR,ILONG
              RETURN
          END IF

          ILG2         = ILONG+2
          FLINV        = 1./FLOAT(ILONG)
          WRKS(ILG2+1) = 0.0

          DO 50 L=1,NL

C             * TRANSFER GR TO WRKS AND DIVIDE BY (2*ILONG).

              DO 20 I=1,ILONG
                  WRKS(I)=GR(I,L)*FLINV
   20         CONTINUE

C             * PERFORM THE FAST FOURIER TRANSFORM.

              N     = ILONG
              ISIGN = -1
              NTYP  =  0

              CALL FOUR2(WRKS,N,1,ISIGN,0)

C             * TRANSFER COMPLEX WAVES 0 TO IR FROM WRKS TO FW.

              NX = -1
              DO 40 J=1,IRP1
                  NX=NX+2
                  FW(J,L)=CMPLX(WRKS(NX  ),WRKS(NX+1  ))
   40         CONTINUE

   50     CONTINUE

      ELSE
C-----------------------------------------------------------------------
C     * CASE 2 - ILONG IS 3 TIMES A POWER OF 2.

C     * MAX VALUE OF IR IS ILONG/3.  WRKS=(ILONG/6+1)*8 COMPLEX WORDS.
C     * THIS SECTION ADDED BY ROGER DALEY, FEB 1976.

          K3             = ILONG/3

          IF (IR.GT.K3) THEN
              WRITE(6,6010) IR,ILONG
              RETURN
          END IF

          FLINV          = 1./FLOAT(ILONG) 
          K6             = K3/2
          KP6            = K6 + 1
          WRKS(14*KP6+1) = 0.0

          CALL PERM(WRKS(4*KP6+1),WRKS(6*KP6+1),AA,BB,KP6,K3)

          DO 150 L=1,NL 

              DO 130 K=1,K3 
                  KP               = (K-1)*3
                  WRKS(8*KP6+K)    = GR(KP+1,L)*FLINV
                  WRKS(10*KP6+K)   = GR(KP+2,L)*FLINV 
                  WRKS(12*KP6+K+1) = GR(KP+3,L)*FLINV
  130        CONTINUE

             WRKS(12*KP6+1)        = WRKS(12*KP6+K3+1)

             N=K3
      ISIGN = -1
      NTYP = 0
      CALL FOUR2(WRKS( 8*KP6+1),N,1,ISIGN,0)
      CALL FOUR2(WRKS(10*KP6+1),N,1,ISIGN,0)
      CALL FOUR2(WRKS(12*KP6+1),N,1,ISIGN,0)

              CALL RCOM (WRKS(4*KP6+1),  WRKS(6*KP6+1), AA, BB,
     1                   WRKS(8*KP6+1),  WRKS(10*KP6+1),
     2                   WRKS(12*KP6+1), WRKS(1), KP6, K3)

              NX = -1
              DO 140 J=1,IRP1
                  NX = NX+2
                  FW(J,L) = CMPLX(WRKS(NX),WRKS(NX+1))
  140         CONTINUE

  150     CONTINUE

      END IF

      RETURN
C-----------------------------------------------------------------------
 6010 FORMAT(' ILLEGAL CALL TO FFWFG2..IR,ILONG=',2I8)
      END 
