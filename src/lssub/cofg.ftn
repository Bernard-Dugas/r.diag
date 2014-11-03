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
C     $Log: cofg.ftn,v $
C     Revision 3.6  2014/09/25 18:31:49  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.5  2013/10/08 00:58:06  bernard
C     Les variables A,P,AZ passent de REAL a REAL(8).
C
C     Revision 3.4  2008/04/28 21:40:20  dugas
C     Corriger l'usage des macros pour r.gppf (passe 2).
C
C     Revision 3.3  1999/11/03 20:45:46  armnrbd
C     Forcer le calcul en Real*8.
C
C     Revision 3.2  1995/11/01 16:56:58  armnrbd
C     De-allouer les appels multi-taches.
C
C     Revision 3.1  1994/11/17  14:12:59  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:55:14  13:55:14  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.1  94/08/09  15:07:02  armnrbd
C     Re-activer les modes de travail hemispheriques.
C     
C     Revision 2.0  93/10/13  13:31:32  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.4  93/03/02  15:31:47  15:31:47  armnrbd (Bernard Dugas)
C     Corriger une erreur survenant lorsque LM=1.
C     
C     Revision 1.3  92/12/21  14:46:50  armnrbd
C     BugFix fonctionnel du niveau 1.2.
C     
C     Revision 1.2  92/12/18  11:44:46  armnrbd
C     Utiliser les symmetries hemispheriques dans les transformees de legendre.
C     
C     Revision 1.1  92/11/13  12:12:53  armnrbd
C     Premier essai: utiliser symmetries hemispheriques.
C     
C     Revision 1.0  92/02/21  11:31:38  armnrbd
C     Initial revision
C     
#     if !defined (nombre_de_taches)
#         define   nombre_de_taches 1
#     endif
      SUBROUTINE cofg ( A,AZ, AP,P,KIND, ILON,ILAT,LA,ILH,
     +                        NP,LM,ILG, CT,LSR,IFAC,W,TRIG)

C *** SUBROUTINE TO TRANSFORM SPECTRAL VALUES TO GRIDPOINT VALUES.

      IMPLICIT none

      INTEGER  ILON,ILAT,LA,ILH,NP,LM,ILG, IFAC(1),LSR(2,1),KIND
      REAL*8   CT(ILG,2,0:NP-1),W(ILG,2,0:NP-1)
      REAL*8   A(2,LA), P(LA,1), AZ(ILON,ILAT)
      REAL*8   AP(2,ILH,ILAT), TRIG(1)

      LOGICAL  HEMI
      INTEGER  I, JP,NLAT,ILATH,LAT0, J1,J2, IUN, M,N
      REAL*8   A1MS,A2MS,A1MN,A2MN, S

CC$    INTEGER  Bidon,mpserv
CC$    EXTERNAL mpserv
      EXTERNAL SFFT
 
      DATA     IUN /   1  /
*--------------------------------------------------------------------
      HEMI  = (KIND.NE.0)
      ILATH = ILAT/2
      IF (HEMI) ILATH = ILAT

***    NLAT IS THE NUMBER OF LATITUDES PAIRS TREATED IN 
***    A PARRALLEL PASS. THIS IS FORCED TO BE EVEN.

      NLAT = (ILATH+NP-1)/NP

CC$doacross local( LAT0, S, I,J1,J2,JP, N,M, A1MS,A2MS,A1MN,A2MN )

      DO 30 JP=0,NP-1

         LAT0 = JP*NLAT+1

         DO 20 J1=LAT0,MIN( ILATH, LAT0+NLAT-1 )

            J2 = ILAT-J1+1

***          LEGENDRE TRANSFORM (SP ==> FR) FOR 
***          LAT( J1 ) AND FOR LAT( J2 ) IN AP.

            DO 10 M=1,LM

               S = -1.0

               A1MS=0.0
               A2MS=0.0

               IF (HEMI)                                       THEN

                  DO N=LSR(1,M+1)-1,LSR(1,M),-1

                     A1MS=A1MS+DBLE( A(1,N) )*P(N,J1)
                     A2MS=A2MS+DBLE( A(2,N) )*P(N,J1)

                  END DO

                  AP(1,M,J1)=A1MS
                  AP(2,M,J1)=A2MS

               ELSE

                  A1MN=0.0
                  A2MN=0.0

***                DETERMINE THE LAST COEFFICIENT'S SYMMETRY.

                  DO N=LSR(1,M),LSR(1,M+1)
                  S = -S
                  END DO

                  DO N=LSR(1,M+1)-1,LSR(1,M),-1

                     S = -S

                     A1MS=A1MS+DBLE( A(1,N) )*P(N,J1)
                     A2MS=A2MS+DBLE( A(2,N) )*P(N,J1)
                     A1MN=A1MN+DBLE( A(1,N) )*P(N,J1)*S
                     A2MN=A2MN+DBLE( A(2,N) )*P(N,J1)*S

                  END DO

                  AP(1,M,J1)=A1MS
                  AP(2,M,J1)=A2MS
                  AP(1,M,J2)=A1MN
                  AP(2,M,J2)=A2MN

               END IF

   10       CONTINUE

   20    CONTINUE

   30 CONTINUE

      NLAT = (ILAT+NP-1)/NP
      NLAT = NLAT+MOD(NLAT,2)

CC$doacross local( LAT0, S, I,J1,J2,JP )

      DO 80 JP=0,NP-1

         LAT0 = JP*NLAT+1

         DO 70 J1=LAT0,MIN( ILAT, LAT0+NLAT-1 ),2

            J2 = J1+1

            DO M=LM+1,ILG+1-LM
               CT(M,1,JP)=0.0
               CT(M,2,JP)=0.0
            END DO

***          (FOURIER COEFFICIENTS) ...
***           ...   LAT( J1 ) + I*LAT( J2 ) IN CT.

            CT(1,1,JP)=AP(1,1,J1)
            CT(1,2,JP)=AP(1,1,J2)
            DO 50 M=2,LM
               CT(      M,1,JP)=+AP(1,M,J1)-AP(2,M,J2)
               CT(      M,2,JP)=+AP(2,M,J1)+AP(1,M,J2)
               CT(ILG+2-M,1,JP)=+AP(1,M,J1)+AP(2,M,J2)
               CT(ILG+2-M,2,JP)=-AP(2,M,J1)+AP(1,M,J2)
   50       CONTINUE

            CALL SFFT( IUN, ILG,1,IFAC, CT(1,1,JP),CT(1,2,JP),
     +                                   W(1,1,JP), W(1,2,JP),TRIG )

            DO 60 I=1,ILG
               AZ(I,J1)=W(I,1,JP)
               AZ(I,J2)=W(I,2,JP)
   60       CONTINUE

            IF (ILG.LT.ILON)                                   THEN
               AZ(ILG+1,J1) = AZ(1,J1)
               AZ(ILG+1,J2) = AZ(1,J2)
            END IF

   70    CONTINUE

   80 CONTINUE

CC$    Bidon = mpserv('BLOCK',Bidon)

      RETURN
      END



