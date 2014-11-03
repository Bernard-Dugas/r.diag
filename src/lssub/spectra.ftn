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
C     $Log: spectra.ftn,v $
C     Revision 3.5  2014/09/25 18:42:04  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.4  2010/02/10 18:19:40  dugas
C     Remplacer les appel a CVMGT par des appels a MERGE
C
C     Revision 3.3  2000/09/25 04:02:12  armnrbd
C     Utiliser COMPLEX*16 dans la routine SPECTRA.
C
C     Revision 3.2  1997/06/02 13:47:48  armnrbd
C     Ajouter la routine QDTFPC.
C
C     Revision 3.1  1994/11/17  14:14:09  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:56:12  13:56:12  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:32:14  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.1  92/04/24  23:13:59  armnrbd
C     Corriger le traitement de CVMGT.
C     
C     Revision 1.0  92/02/21  11:34:39  armnrbd
C     Initial revision
C     

      SUBROUTINE spectra (F,LM,LR,KTR,LA,LSR,SUM,KD)

***    AUG 08/90 - G.J.BOER, F.MAJAESS(ADD LR, KTR AND RESTRUCTURE DO LOOPS)
***    FEB XX/80 - S. LAMBERT

***     SPECTRA READS IN A COMPLEX-VALUED SPECTRAL FILE F(LA) AND RETURNS
***      THE SUM OVER N/M AS A FUNCTION OF M/N IF (KD .EQ./.NE. 2).
***      NOTE:  KTR=0/2 INDICATES RHOMBOIDAL/TRIANGULAR TRUNCATION.

      IMPLICIT none

      INTEGER    LR,LM,LA,KTR,KD,LSR(2,1),
     +           K,M,MF,ML,N,NF,NL,NMAX
      COMPLEX*16 F(LA),SUM(1)

*-----------------------------------------------------------------------
      IF (LR.LT.2 .OR.
     +    LM.LT.2)    CALL                         XIT(' Spectra',-1 )
      IF (LR.LT.LM)   CALL                         XIT(' Spectra',-2 )

***    CHECK IF THE SUM IS OVER M OR N.

      IF (KD.EQ.2)                                             THEN
 
***        SUM OVER N IS REQUESTED.
 
          DO 200 M=1,LM
              SUM(M) = 0.0
              NF     = M
              NL     = MERGE( LR+M-1, LR, KTR.EQ.0 )
              DO 100 N=NF,NL
                  K      = LSR(1,M)+(N-M)
                  SUM(M) = SUM(M)+F(K)
 100          CONTINUE
 200      CONTINUE
 
      ELSE
 
***        SUM OVER M IS REQUESTED.
 
          NMAX = MERGE( LM+LR-1, LR, KTR.EQ.0 )
 
          DO 400 N=1,NMAX
              SUM(N) = 0.0
              MF     = MERGE( MAX0( 1,(N-LR+1) ), 1, KTR.EQ.0 )
              ML     = MIN0( N,LM )
              DO 300 M=MF,ML
                  K      = LSR(1,M)+(N-M)
                  SUM(N) = SUM(N)+F(K)
 300          CONTINUE
 400      CONTINUE
 
      END IF
 
      RETURN
*-------------------------------------------------------------------

      END
      SUBROUTINE qdtfpc (Q,D,P,C,LA,LSR,LM,ILEV,KASE)

***    APR 11/79 - J.D.HENDERSON 

***    Q(LA,ILEV) = SPECTRAL VORTICITY 
***    D(LA,ILEV) = SPECTRAL DIVERGENCE
***    P(LA,ILEV) = SPECTRAL STREAMFUNCTION
***    C(LA,ILEV) = SPECTRAL VELOCITY POTENTIAL

***    KASE = +1 CONVERTS P,C TO Q,D.
***    KASE = -1 CONVERTS Q,D TO P,C.

      IMPLICIT none

      INTEGER  LA,LSR(2,1),LM,ILEV,KASE
      COMPLEX  Q(LA,ILEV),D(LA,ILEV),P(LA,ILEV),C(LA,ILEV) 

***    LOCAL VARIABLES.

      INTEGER  K,L,M,KL,KR,NS
      REAL     FNS1,RFNS1

*-------------------------------------------------------------------
      IF (KASE.GE.0)                                           THEN

***        KASE = +1 CONVERTS P,C TO Q,D.

          DO  100 L=1,ILEV 
              DO  M=1,LM 
                  KL = LSR(1,M) 
                  KR = LSR(1,M+1)-1 
                  DO  K=KL,KR
                      NS     = (M-1)+(K-KL) 
                      FNS1   = FLOAT(NS*(NS+1)) 
                      Q(K,L) =-FNS1*P(K,L) 
                      D(K,L) =-FNS1*C(K,L) 
                  END DO
              END DO
  100     CONTINUE
          
      ELSE

***        KASE = -1 CONVERTS Q,D TO P,C.

          RFNS1=0.
          DO  200 L=1,ILEV 
              DO  M=1,LM 
                  KL = LSR(1,M) 
                  KR = LSR(1,M+1)-1 
                  DO  K=KL,KR
                      NS=(M-1)+(K-KL) 
                      IF(NS.GT.0) RFNS1=1./FLOAT(NS*(NS+1)) 
                      P(K,L)=-RFNS1*Q(K,L)
                      C(K,L)=-RFNS1*D(K,L)
                  END DO
              END DO
  200     CONTINUE

      END IF

      RETURN
*-------------------------------------------------------------------

      END
