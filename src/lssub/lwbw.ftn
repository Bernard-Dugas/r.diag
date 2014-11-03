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
C     $Log: lwbw.ftn,v $
C     Revision 3.6  2014/09/25 18:42:03  dugas
C     Inserer le texte de la licence LPGPL 2.1 pour R.DIAG
C
C     Revision 3.5  2013/10/08 01:12:51  bernard
C      - La variables W est maintenant declaree en REAL(8)
C
C     Revision 3.4  1998/04/07 20:51:32  armnrbd
C     Corriger la copie de la derniere latitude si NLON1 est impair.
C
C     Revision 3.3  1997/06/02  14:27:58  armnrbd
C     Ajouter la routine IMAVRAI.
C
C     Revision 3.2  1996/02/06  18:15:49  armnrbd
C     Introduire NLON1 et definir NLON selon la parite de NLON1.
C
C     Revision 3.1  1994/11/17  14:13:44  armnrbd
C     Messages informatifs quand au passage de la version 2.x a 3.1...
C     1) Les espaces en debut des noms de variables de sont plus pertinents.
C     2) Les grilles complexes de type CMPL sont maintenant supportees.
C     3) Les fichiers SQI sont reconnus, lus et ecrit directements.
C     4) Plusieurs nouvelles cles sont disponibles au demarrage.
C
C     Revision 3.0  94/11/17  13:55:51  13:55:51  armnrbd (Bernard Dugas)
C     *** empty log message ***
C     
C     Revision 2.0  93/10/13  13:31:58  armnrbd
C     Premiere version compatible HP-UX.
C     
C     Revision 1.0  92/02/21  11:33:48  armnrbd
C     Initial revision
C     

      SUBROUTINE lwbw (W,NLON1,NLAT,NLEV,CL,KIND)

C     * NOV 20/90 - B.DUGAS, RPN. (PORT TO CY920/ IRIX F77)
C     * FEB 26/80 - J.D.HENDERSON 

C     * CONVERSION BETWEEN REAL WIND COMPONENTS (LITTLE WINDS) AND
C     *  MODEL WIND COMPONENTS (BIG WINDS). 

C     * W    CONTAINS WIND COMPONENT TO BE CONVERTED.
C     * CL   CONTAINS COS(LAT) S TO N.

C     * THE TYPE OF CONVERSION DEPENDS ON THE VALUE OF KIND.
C     * +1 = LITTLE TO BIG,  -1 = BIG TO LITTLE,
C     * OTHERWISE RETURN WITH NO CONVERSION.

      IMPLICIT NONE

      INTEGER KIND,NLON1,NLAT,NLEV, I,J,L, NLON
      REAL*8  W(NLON1,NLAT,1)
      REAL*8  CL(1), A,CON

C---------------------------------------------------------------------
      A = 6.37122D6

                               NLON = NLON1-1
      IF (MOD( NLON1,2 ).EQ.0) NLON = NLON1

      IF (NLON1.NE.1)                                          THEN

          IF (KIND.EQ.+1)                                      THEN

C             * KIND=1 CONVERTS LITTLE WINDS TO BIG WINDS
C             * BY MULTIPLYING BY COS(LAT)/(EARTH RADIUS).

              DO 210 L=1,NLEV 
                  DO 210 J=1,NLAT 

                      CON = CL(J)/A 
                      DO 209 I=1,NLON
                          W(I,J,L) = W(I,J,L)*CON 
  209                 CONTINUE

                      IF (NLON.NE.NLON1)
     +                W(NLON1,J,L)=W(1,J,L)

  210         CONTINUE

          ELSE IF (KIND.EQ.-1)                                 THEN

C             * KIND=-1 CONVERTS BIG WINDS TO LITTLE WINDS
C             * BY MULTIPLYING BY EARTH RADIUS/COS(LAT).

              DO 310 L=1,NLEV 
                  DO 310 J=1,NLAT 

                      CON = A/CL(J) 
                      DO 309 I=1,NLON
                          W(I,J,L) = W(I,J,L)*CON 
  309                 CONTINUE

                      IF (NLON.NE.NLON1)   
     +                W(NLON1,J,L) = W(1,J,L)

  310         CONTINUE

          END IF

      ELSE

C         * SPECIAL CASE WHERE THE FIELDS ARE ZONAL CROSS-SECTIONS.

          IF (KIND.EQ.+1)                                      THEN

C             * KIND=1 CONVERTS LITTLE WINDS TO BIG WINDS
C             * BY MULTIPLYING BY COS(LAT)/(EARTH RADIUS).

              DO 410 L=1,NLEV 
                  DO 410 J=1,NLAT 

                      W(1,J,L) = CL(J)*W(1,J,L)/A

  410         CONTINUE

          ELSE IF (KIND.EQ.-1)                                 THEN

C             * KIND=-1 CONVERTS BIG WINDS TO LITTLE WINDS
C             * BY MULTIPLYING BY EARTH RADIUS/COS(LAT).

              DO 510 L=1,NLEV 
                  DO 510 J=1,NLAT 

                      W(1,J,L) = A*W(1,J,L)/CL(J)

  510         CONTINUE

          END IF
    
      END IF

      RETURN
*---------------------------------------------------------------- 

      END 
      SUBROUTINE imavrai (UG,VG,PSDLG,PSDPG, COSJ,ILG,LON,ILEV,A) 

***    JUL 14/92 - E. CHAN (ADD REAL*8 DECLARATIONS)
***    JAN 12/88 - R.LAPRISE.

***    CONVERT WIND IMAGES AND D(PS)/DX AND /DY TO REAL VALUES.

      IMPLICIT none

      INTEGER  ILG,LON,ILEV, I,L
      REAL     UG(ILG,ILEV),VG(ILG,ILEV),PSDLG(ILG),PSDPG(ILG),A
      REAL*8   COSJ

*---------------------------------------------------------------- 
      DO 100 L=1,ILEV 
         DO I=1,LON 
            UG(I,L)   = UG(I,L)     *(A/COSJ) 
            VG(I,L)   = VG(I,L)     *(A/COSJ) 
         END DO
  100 CONTINUE

      DO 200 I=1,LON
            PSDLG(I)  = PSDLG(I)*(1./(A*COSJ))
            PSDPG(I)  = PSDPG(I)*(1./(A*COSJ))
  200 CONTINUE

      RETURN
*-----------------------------------------------------------------------

      END
