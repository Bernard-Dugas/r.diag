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
C     $Log: ffgfw2.ftn,v $
C     Revision 3.2  2014/09/25 18:42:02  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.1  1994/11/17 14:13:12  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:55:25  13:55:25  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:31:39  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.0  92/02/21  11:32:30  armnrbd
C     Initial revision
C     

      SUBROUTINE ffgfw2 (GR,LGR,FW,LFW,IR,ILONG,WRKS,NL)

C     *****   FEB 1976  -  JOHN D. HENDERSON  ****

C     * NOV 24/82 - R.LAPRISE, S.J.LAMBERT.

C     * USE THE FORTRAN FAST FOURIER TRANSFORM TO PERFORM A WAVE TO 
C     * GRID FOURIER TRANSFORM OF MULTI-LEVEL FIELD.

C     * NOTE... GR AND FW ARE GENERALLY EQUIVALENCED IN CALLING ROUTINE

C     * PARAMETERS...
C     * GR    = MULTI-LEVEL ARRAY TO RETURN GRID POINT FIELDS. 
C     * LGR   = FIRST DIMENSION OF GR.
C     * FW    = MULTI-LEVEL COMPLEX ARRAY CARRYING FOURIER AMPLITUDES.
C     * LFW   = FIRST DIMENSION OF FW.
C     * IR    = MAXIMUM E.W. WAVENUMBER TO BE USED IN THE TRANSFORM.
C     * ILONG = NUMBER OF LONGITUDES ON THE TRANSFORM GRID.
C     * WRKS  = WORK FIELD.
C     * NL    = NUMBER OF VERTICAL LEVELS.

      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)

      COMPLEX    FW(LFW,NL),WRKS(1)
      DIMENSION  GR(LGR,NL)
      COMPLEX    AA,BB,C

C-----------------------------------------------------------------------
      IRP1 = IR+1

      IF (MOD(ILONG,3).NE.0) THEN
C-----------------------------------------------------------------------
C     * CASE 1 - ILONG IS A POWER OF 2. 

C     * MAX IR IS ILONG/2.  WRKS=(ILONG+1) COMPLEX WORDS.

          IF (IR.GT.ILONG/2) THEN
              WRITE(6,6010) IR,ILONG
              RETURN
          END IF

          ILG2          = ILONG+2
          NWMAX         = ILONG/2  + 1
          WRKS(NWMAX+1) = (0.,0.)

          DO 50 L=1,NL

C         * COPY WAVES 0 TO IR TO WRKS. FILL TO WAVE ILONG/2 WITH ZEROS.

              DO 20 J=IRP1,NWMAX
                  WRKS(J)=(0.,0.)
   20         CONTINUE

              DO 30 J=1,IRP1
                  WRKS(J)=FW(J,L)
   30         CONTINUE

C             * PERFORM THE FAST FOURIER TRANSFORM.

              N     = ILONG
              ISIGN =  1
              NTYP  = -1

              CALL FOUR2(WRKS,N,1,ISIGN,-1)

C             * COPY ILONG GRID POINTS FROM WRKS TO GR.

              NX=0
              DO 40 I=1,ILONG,2
                  NX=NX+1
                  GR(I  ,L)= REAL(WRKS(NX))
                  GR(I+1,L)=AIMAG(WRKS(NX))
   40         CONTINUE

   50     CONTINUE

      ELSE
C-----------------------------------------------------------------------
C     * CASE 2 - ILONG IS 3 TIMES A POWER OF 2.

C     * MAX VALUE OF IR IS ILONG/3.  WRKS=(ILONG/6+1)*8 COMPLEX WORDS.
C     * THIS SECTION ADDED BY ROGER DALEY, FEB 1976.

          K3            = ILONG/3

          IF (IR.GT.K3) THEN
              WRITE(6,6010) IR,ILONG
              RETURN
          END IF

          K6            = K3/2
          KP6           = K6 + 1
          NWMAX         = 2*KP6 

          WRKS(7*KP6+1) = (0.,0.) 

          CALL PERM(WRKS(2*KP6+1),WRKS(3*KP6+1),AA,BB,KP6,K3)

          DO 160 L=1,NL 

              DO 120 J=IRP1,NWMAX
                  WRKS(J) = (0.,0.)
  120         CONTINUE

              DO 130 J=1,IRP1
                  WRKS(J) = FW(J,L)
  130         CONTINUE

              DO 140 K=1,KP6
                  CC            = CONJG(WRKS(K3+2-K))
                  WRKS(4*KP6+K) =  WRKS(K) + C
                  WRKS(5*KP6+K) = (WRKS(K) + BB*CC)*WRKS(2*KP6+K)
                  WRKS(6*KP6+K) = (WRKS(K) + AA*CC)*WRKS(3*KP6+K)
  140         CONTINUE

              N     = K3
              ISIGN =  1
              NTYP  = -1

              CALL FOUR2(WRKS(4*KP6+1),N,1,ISIGN,-1)
              CALL FOUR2(WRKS(5*KP6+1),N,1,ISIGN,-1)
              CALL FOUR2(WRKS(6*KP6+1),N,1,ISIGN,-1)

              WRKS(7*KP6) = WRKS(6*KP6+1)

              DO 150 K=1,K6 
                  KP         = (K-1)*6
                  GR(KP+1,L) =  REAL(WRKS(4*KP6+K))
                  GR(KP+2,L) =  REAL(WRKS(5*KP6+K))
                  GR(KP+3,L) = AIMAG(WRKS(6*KP6+K))
                  GR(KP+4,L) = AIMAG(WRKS(4*KP6+K))
                  GR(KP+5,L) = AIMAG(WRKS(5*KP6+K))
                  GR(KP+6,L) =  REAL(WRKS(6*KP6+K+1)) 
  150         CONTINUE

  160     CONTINUE

      END IF

      RETURN
C-----------------------------------------------------------------------
 6010 FORMAT(' ILLEGAL CALL TO FFGFW2..IR,ILONG=',2I8)
      END 
