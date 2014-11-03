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
C     $Log: gcof.ftn,v $
C     Revision 3.6  2014/10/16 12:00:39  dugas
C     Modifications commandees par la compilation avec GFORTRAN et IFORT.
C
C     Revision 3.5  2014/09/25 18:42:02  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.4  2013/10/08 00:58:06  bernard
C     Les variables A,P,AZ passent de REAL a REAL(8).
C
C     Revision 3.3  2008/04/28 21:45:28  dugas
C     Corriger l'usage des macros pour r.gppf (passe 2).
C
C     Revision 3.2  1995/11/01 20:04:03  armnrbd
C     De-allouer les appels multi-taches.
C
C     Revision 3.1  1994/11/17  14:13:27  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:55:34  13:55:34  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.3  94/10/12  14:08:37  14:08:37  armnrbd (Bernard Dugas)
C     Corriger la mise a zero dans le cas anti-symmetrique hemispherique.
C     
C     Revision 2.2  94/08/24  12:10:50  armnrbd
C     Mettre a zero les composantes symmetriques ou anti-symmetriques
C     du nombre zonal 0 dans le cas de transformees de parites opposees.
C     
C     Revision 2.1  94/08/09  15:07:08  armnrbd
C     Re-activer les modes de travail hemispheriques.
C     
C     Revision 2.0  93/10/13  13:31:45  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.4  92/12/21  15:10:13  15:10:13  armnrbd (Bernard Dugas)
C     BugFix fonctionnel du niveau 1.3.
C     
C     Revision 1.3  92/12/17  22:09:52  armnrbd
C     Utiliser les symmetries hemispheriques dans les transformees de legendre.
C     
C     Revision 1.2  92/11/15  22:34:00  armnrbd
C     Correction a la section zonale.
C     
C     Revision 1.1  92/11/13  12:15:19  armnrbd
C     Premier essai: utiliser symmetries hemispheriques.
C     
C     Revision 1.0  92/02/21  11:33:07  armnrbd
C     Initial revision
C     
#     if !defined (nombre_de_taches)
#         define   nombre_de_taches 1
#     endif
      SUBROUTINE gcof ( AZ,A, AP,P,KIND,    ILON,ILAT,LA,ILH,
     +                        NP,LM,ILG,CT, LSR,IFAC,W,TRIG,GUSW)

***    SUBROUTINE TO TRANSFORM GRIDPOINT VALUES TO SPECTRAL VALUES.

      IMPLICIT none

      INTEGER  ILON,ILAT,LA,ILH,NP,LM,ILG,IFAC(1),LSR(2,*),KIND
      REAL*8   A(2,LA),AZ(ILON,ILAT),P(LA,ILAT)
      REAL*8   CT(ILG,2,0:NP-1),W(ILG,2,0:NP-1),AP(2,ILH,ILAT)
      REAL*8   GUSW(ILAT),TRIG(*)

      LOGICAL  HEMI
      INTEGER  I, JP, M,N, LAT0,NLAT,ILATH
      INTEGER  MF,MI,NLD,NLF,NLI,N0,J1,J2, MUN
      REAL*8   A1MS,A2MS,A1MN,A2MN,PN1,S,OVIM

CC$    INTEGER  Bidon,mpserv
CC$    EXTERNAL mpserv
      EXTERNAL SFFT

      DATA     MUN / -1 /
*--------------------------------------------------------------------
      HEMI  = (KIND.NE.0)
      ILATH = ILAT/2
      IF (HEMI) ILATH = ILAT

      IF (ILON.NE.1)                                           THEN

         OVIM = 1.0/DBLE( 2*ILG )

***       NLAT IS THE NUMBER OF LATITUDES
***       TREATED IN A PARRALLEL PASS.

         NLAT = (ILAT+NP-1)/NP
         NLAT = NLAT+MOD(NLAT,2)

CC$mp_schedtype=INTERLEAVE
CC$doacross local( I, JP, N,M, J1,J2, LAT0 )

         DO 40 JP=0,NP-1

            LAT0 = JP*NLAT+1

            DO 30 J1=LAT0,MIN( ILAT, LAT0+NLAT-1 ),2

               J2 = J1+1

***             FOURIER TRANSFORM (GG ==> FR) FOR 
***             LAT( J1 ) AND FOR LAT( J2 ) IN CT.

               DO 10 I=1,ILG
                  CT(I,1,JP)=AZ(I,J1)*OVIM
                  CT(I,2,JP)=AZ(I,J2)*OVIM
   10          CONTINUE

               CALL SFFT( MUN, ILG,1,IFAC, CT(1,1,JP),CT(1,2,JP),
     +                                      W(1,1,JP), W(1,2,JP),TRIG )

               AP(1,1,J1) = 2.0*W(1,1,JP)
               AP(1,1,J2) = 2.0*W(1,2,JP)

               AP(2,1,J1) = 0.0
               AP(2,1,J2) = 0.0

               DO 20 M=2,LM
                  AP(1,M,J1)=+W(M,1,JP)+W(ILG+2-M,1,JP)
                  AP(2,M,J1)=+W(M,2,JP)-W(ILG+2-M,2,JP)
                  AP(1,M,J2)=+W(M,2,JP)+W(ILG+2-M,2,JP)
                  AP(2,M,J2)=-W(M,1,JP)+W(ILG+2-M,1,JP)
   20          CONTINUE

   30       CONTINUE

   40    CONTINUE


         DO 70 MI=1,LM,NP
         MF = MIN( MI+NP-1,LM )

CC$doacross local(M,N,N0,J1,J2,NLD,NLI,NLF,A1MS,A2MS,A1MN,A2MN,S,PN1)

         DO 70 M=MI,MF
            N0 = M-MI

            NLI = LSR(1,M)
            NLF = LSR(1,M+1)-1
            NLD = NLF-NLI+1

            DO N=1,NLD
               W(N,1,N0) = 0.
               W(N,2,N0) = 0.
            END DO

            DO 60 J1=1,ILATH

               J2  = ILAT-J1+1
               S = -1.0D0

               A1MS = AP(1,M,J1)*GUSW(J1)
               A2MS = AP(2,M,J1)*GUSW(J1)

               IF (HEMI)                                       THEN

                  DO N=1,NLD

                     PN1 =  P(N+NLI-1,J1)

                     W(N,1,N0)=W(N,1,N0)+ A1MS*PN1
                     W(N,2,N0)=W(N,2,N0)+ A2MS*PN1

                  END DO

               ELSE

                  A1MN = AP(1,M,J2)*GUSW(J2)
                  A2MN = AP(2,M,J2)*GUSW(J2)

                  DO N=1,NLD

                     S = -S
                     PN1 =  P(N+NLI-1,J1)

                     W(N,1,N0)=W(N,1,N0)+(A1MS+S*A1MN)*PN1
                     W(N,2,N0)=W(N,2,N0)+(A2MS+S*A2MN)*PN1

                  END DO

               END IF

   60       CONTINUE

            DO N=1,NLD
               A(1,N+NLI-1) = W(N,1,N0)
               A(2,N+NLI-1) = W(N,2,N0)
            END DO

   70    CONTINUE

CC$       Bidon = mpserv('BLOCK',Bidon)

      ELSE

***       SPECIAL CASE WHERE DATA CONSIST OF ZONAL CROSS-SECTIONS.

         S = -1.0

         DO 1060 N=1,LSR(1,2)-1

            S    = -S
            A1MS = 0.D0

            IF (HEMI)                                          THEN

               DO J1=1,ILATH
                  A1MS = A1MS+AZ(1,J1)*P(N,J1)*GUSW(J1)
               END DO

            ELSE

               DO J1=1,ILATH
                  J2   = ILAT-J1+1
                  A1MS = A1MS+AZ(1,J1)*P(N,J1)*GUSW(J1)
     +                 +    S*AZ(1,J2)*P(N,J1)*GUSW(J2)
               END DO

            END IF

            A(1,N)=A1MS
            A(2,N)=0.0

 1060    CONTINUE

      END IF

      IF (KIND.GT.0)                                           THEN

***        ZERO-OUT THE ANTI-SYMMETRIC ZONAL COMPONENT.

          DO  N=2,LSR(1,2)-1,2
              A(1,N) = 0.0
          END DO

      ELSE IF (KIND.LT.0)                                      THEN

***        ZERO-OUT THE SYMMETRIC ZONAL COMPONENT.

          DO  N=1,LSR(1,2)-1,2
              A(1,N) = 0.0
          END DO

      END IF

      RETURN
      END
