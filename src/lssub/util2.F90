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
!
! Description: Ceci est le trop plein de util.ftn en F90
!
!     Log: util2.F90
!     Revision 1.3  2016/10/25 14:44  dugas
!     - Ajouter le texte de la licence LGPL (oubli inexcusable)
!     - Ajouter la routine DIAG_CCARD (clone de CCARD en F90)
!
!     Revision 1.1+1.2  2014-12-09 14:14:06 dugas
!     - Edition d'un commentaire.
!     - Simplifier les fonctions INTEFLT et ME32O64.

!     Revision 1.0  2014/12/03 23:24:37  dugas
!     Version initiale
!

subroutine messys( string )
#  if defined (__HOS_AIX__)
   USE XLFUTILITY_EXTNAME
   character * (*) string
   character(len=512) ligne
   integer ier
   ier = ierrno( )
   write(ligne,'("SYSERROR no. ",I5)') ier
#  else
#  if defined (__INTEL_COMPILER_UPDATE)
   USE IFCORE, only: GERROR
#  endif
   character * (*) string
   character(len=512) ligne
   call gerror( ligne )
#  endif
   write(6,'(A)') string
   write(6,'(A)') ligne

   return

end subroutine messys

function CHKENDI() result(outval)

   use ISO_C_BINDING  

   implicit none

   ! Determines the endian-ess of a computer
   ! Version re-coded on Nov 6, 2014 (after M. Valin)

   integer :: outval

   type(C_PTR) :: ptr

   integer (kind=2), dimension(:), pointer :: I2
   integer (kind=4),                target :: I4

   I4 = 1 ; ptr = c_loc(I4) ; call c_f_pointer( ptr,I2,[2] )

   outval = 1 - I2(1) ! = 0 (Little) or 1 (Big) endian
   ! if ( I2(1) == 1 ) outval = 0 ! Little endian machine
   ! if ( I2(1) == 0 ) outval = 1 ! Big endian machine

   return

end function CHKENDI

function ME32O64() result(outval)

   implicit none

   integer :: outval
 
   ! Tell calling program whether the current environment is
   ! configured with default 32-bit integers (2) or with
   ! default 64-bit integers (1)

   ! Original version by A.J. Stacey -MARCH 23,1992-
   ! Re-coded (following M. Valin) by B. Dugas, Nov 2014

   integer, dimension(2) :: INUM

   outval = -1

   if (loc( INUM(2) ) - loc( INUM(1) ) == 4) then
      outval = 2 ! 32-bits integers
   else if (loc( INUM(2) ) - loc( INUM(1) ) == 8) then
      outval = 1 ! 64-bits integers
   end if
 
   return

end function ME32O64

function INTEFLT() result(outval)

   implicit none

   integer :: outval
 
   ! This routine checks if integer word size is the same
   ! as the real word size. This is subtly different from
   ! the output of "ME32O64".

   ! Original version by A.J. Stacey - March 1992
   ! Re-coded (following M. Valin) by B. Dugas, Nov 2014

   real,    dimension(2) :: ANUM
   integer, dimension(2) :: INUM

   if ((loc( INUM(2) ) - loc( INUM(1) )) == &
       (loc( ANUM(2) ) - loc( ANUM(1) ))) then
      outval = 1 ! default integer and real sizes are identical
   else
      outval = 2 ! default integer and real sizes are different
   end if
 
   return

end function INTEFLT


#if !defined (NKLEMAX)
#define NKLEMAX 1024
#endif

!**S/P DIAG_CCARD - RECUPERATION DES PARAMETRES D'APPEL A UN PROGRAMME

SUBROUTINE DIAG_CCARD( INCLE,DEF,VAL,N,II )

   IMPLICIT NONE

   INTEGER N,II
   CHARACTER * (*) INCLE(N), DEF(N), VAL(N)

!ADAPTATION: JAMES CAVEEN (1991)
!            B. Dugas (Version F90, octobre 2016)
!
!OBJET(DIAG_CCARD)- RECUPERATION DES PARAMETRES D'APPEL A UN PROGRAMME
!                   SI LE NOM D'UNE CLEF APPARAIT SEUL LORS DE LA SEQUENCE
!                   D'APPEL, ON LUI ATTRIBUE LA VALEUR CONTENUE DANS DEF
!                   SI LE NOM APPARAIT, SUIVI D'UNE OU DES VALEURS, ON LUI
!                   ATTRIBUE CES VALEURS.  SINON, LA CLEF PREND POUR VALEUR
!                   LE CONTENU INITIAL DE VAL.
!
!                   ON PEUT DONNER UNE VALEUR A UNE CLEF DES FACONS
!                   SUIVANTES:
!
!                   -CLEF VALEUR
!                   -CLEF VAL1:VAL2:VAL3:...
!                   -CLEF 2:3:-4
!                   -CLEF =-2:3:4
!                   -CLEF =-22
!                   -CLE=-22:33:...           
!      
!                   LORS DE LA SEQUENCE D'APPEL, TOUS LES PARAMETRES
!                   PRECEDANT LA PREMIERE CLEF RENCONTREE ET
!                   SUIVANT LES SIGNES -- SONT TRAITES EN MODE POSITIONEL.
!
!                   EXEMPLES D'UTILISATION:
!
!                   PROGNOM -CLEF1 -CLEF2 VAL2 -- POS1 POS2
!                   PROGNOM POS1 POS2 -CLEF1 -- POS3 POS4 ...
!                   PROGNOM -CLEF1 =-12:33 -CLEF2 12 =-34
!
!                   LORSQUE LE PREMIER PARAMETRE PASSE = -H/-h, CCARD IMPRIME
!                   SUR L'UNITE 6 LA SEQUENCE D'APPEL AU PROGRAMME.
!                 
!                   LA SEQUENCE D'APPEL ET LES VALEURS COURANTES DES CLEFS
!                   SERONT IMPRIMEES S'IL Y A ERREUR LORS DE L'APPEL.
!
!
!ARGUMENTS:
!       INCLE     ENTREE     -  NOM DES DIFFERENTES CLEFS A INITIALISER
!       DEF          "       -  DEUXIEME VALEUR DE DEFAUT
!       VAL       SORTIE     -  VALEUR FINALE ATTRIBUEE AUX CLEFS
!       N         ENTREE     -  NOMBRE DE CLEFS A INITIALISER
!       II        SORTIE     -  NOMBRE DE PARAMETRES POSITIONNELS
!                               (I.E. ASSOCIES AUX CLEFS "-")
!                            -  SI II > 0, ON RETOURNE LE NOMBRE DE 
!                               PARAMETRES POSITIONELS.
!                            -  SI II = -111 ou +111, ON ARRETE DES LA PREMIERE
!                               ERREUR
!                            -  SI II <=0 & <> -111, ON NE RETOURNE PAS
!                               LE NOMBRE D'ARGUMENTS POSITIONELS ET ON
!                               CONTINUE MALGRE LES ERREURS
!*
      
   INTEGER, external :: LONGUEUR
   EXTERNAL LOW2UP, XIT

   INTEGER I,POS,J,POSC,IINDEX,POSMOIN,POSMOINC

   ! TYPE DE CLE -1=MINUS, 1=MAJUS, 0=PAREIL
   INTEGER CLTYPE(NKLEMAX) 

   CHARACTER(LEN=50)   CLEUP, CLE(NKLEMAX), cletemp
   CHARACTER(LEN=8192) ARGUP, ARG, CCARD_ARGS, ARGtemp
   CHARACTER(LEN=60)   Keyname
   LOGICAL PASFINI, PLANTE
   INTEGER STATUS, INTERNE
   Integer L_argenv

   ! SI II = -111, LE PROGRAMME ARRETE DES LA PREMIERE 
   ! ERREUR RENCONTREE SINON, ON CONTINUE

   plante = (abs(ii) .eq. 111)

   call get_environment_variable('CCARD_OPT',ccard_args,L_argenv)

   IF (L_argenv .gt. 0) THEN
      IF (ccard_args(1:L_argenv) .eq. 'ABORT') THEN
         plante = .true.
      ENDIF
   ENDIF

   call get_environment_variable('CCARD_ARGS',ccard_args,L_argenv)

   ! INITIALISER LE VECTEUR DE TYPE DE CLEFS

   DO I = 1, N

      cletemp = incle(i)
      IF (TRIM(incle(i)) .ne. TRIM(cletemp)) THEN
         write(6,777) 'CCARD erreur: nom de la cle #',i,' > limite de 50 caracteres'
         if (plante) call xit('  Ccard ',-21)
      ENDIF
      CALL LOW2UP(INCLE(I),CLEUP)
      DO J =50,1,-1
         IF(CLEUP(J:J) .NE. ' ') THEN
            IINDEX = J
            EXIT
         ENDIF
      ENDDO

      IF(CLEUP(IINDEX:IINDEX) .EQ. '.') THEN
         CLTYPE(I) = 0
         CLEUP(IINDEX:IINDEX) = ' '
      ELSE IF (CLEUP(IINDEX:IINDEX) .EQ. '_') THEN
         CLTYPE(I) = -1
         CLEUP(IINDEX:IINDEX) = ' '
      ELSE IF(CLEUP(IINDEX:IINDEX) .EQ. ':') THEN
         CLTYPE(I) = 2
         CLEUP(IINDEX:IINDEX) = ' '
      ELSE
         CLTYPE(I) = 1
      ENDIF
      CLE(I) = CLEUP

   ENDDO

   ! TROUVER LA POSITION DE LA PREMIERE CLEF - DANS LA LISTE CLE

   POSMOIN = 0 ; POSMOINC = 0

   DO J=1,N
      IF (CLE(J) .EQ. '-') THEN
         POSMOINC = J
         POSMOIN  = J
         EXIT
      ENDIF
   ENDDO

   POS  = 0
   POSC = 0
   INTERNE   = 0
   PASFINI = .TRUE.

   DO WHILE (PASFINI)

      if (L_argenv .gt. 0) then
         status = qqqoenv(arg,ccard_args,L_argenv)
      else
         STATUS = QQQOBM(ARG,cltype(pos))
      endif

      IF(STATUS .EQ. 1) THEN
         CALL LOW2UP(ARG,ARGUP)
         DO J = 1,N
            CLEUP = CLE(J)
            IF (trim(ARGUP) .EQ. trim(CLEUP)) THEN
               POSC = J
               POS = J
               CALL QQQTRNS(VAL(POS), DEF(POS),CLTYPE(POS))
               ! print *,'Debug+ ccard cleup=',cleup,' val(pos)=',val(pos)
               write(keyname,109) '%%'//trim(cle(posc)),posc-pos,'%%'
               call c_set_appl_var(Keyname,val(pos))
               EXIT
            ENDIF
         ENDDO
         IF (J == N+1) THEN
            PRINT *,' *** ERREUR: cle ',TRIM( ARGUP ),' non reconnue'
            CALL QQQSAP(CLE,DEF,VAL,N)
            IF(PLANTE) CALL XIT('  Ccard ',-22 )
         ENDIF
      ELSE IF(STATUS .EQ. 2) THEN
         IF (POSC .NE. 0 .AND. CLE(POS) .EQ. CLE(POSC)) THEN
            CALL QQQTRNS(ARGtemp,ARG,CLTYPE(POSC))
            VAL(POSC) = ARGtemp
            ! print *,'Debug+ ccard cle(posc)=',cle(posc),' ARGtemp=',trim(ARGtemp)
            ! print *,'Debug+ ccard val(posc)=',val(posc)
            write(keyname,109) '%%'//trim(cle(posc)),posc-pos,'%%'
            call c_set_appl_var(Keyname,ARGtemp)
            POSC = POSC + 1
         ELSE
            PRINT *,'ARG=',TRIM(ARG)
            WRITE(6,'("POS,POSC=",2I4)') POS,POSC
            PRINT *,' *** ERREUR: debordement de liste (valeur argument)'
            CALL QQQSAP(CLE,DEF,VAL,N)
            IF(PLANTE) CALL XIT('  CCARD ',-23)
         ENDIF
      ELSE IF(STATUS .EQ. 3) THEN
         IF (POSMOINC .NE. 0 .AND. CLE(POSMOIN) .EQ. CLE(POSMOINC)) THEN
            CALL QQQTRNS(ARGtemp,ARG,CLTYPE(POSMOINC))
            VAL(POSMOINC) = ARGtemp
            write(keyname,109) '%%'//trim(cle(posmoinc)),posmoinc-pos,'%%'
            call c_set_appl_var(Keyname,ARGtemp)
         ELSE
            PRINT *,'ARG=',TRIM(ARG)
            WRITE(6,'("POS,POSC,POSMOIN,POSMOINC=",4I4)') POS,POSC,POSMOIN,POSMOINC
            IF (POSMOIN > 0 .AND. POSMOINC > 0) &
                WRITE(6,'("CLE(POSMOIN),CLE(POSMOINC)=",2A)') TRIM(CLE(POSMOIN)), TRIM(CLE(POSMOINC))
            PRINT *,' *** ERREUR: debordement de liste'
            PRINT *,'         ou  mode positionnel non permis'
            CALL QQQSAP(CLE,DEF,VAL,N)
            IF(PLANTE) CALL XIT('  CCARD ',-24)
         ENDIF
         POSMOINC = POSMOINC + 1
         INTERNE = INTERNE + 1
      ELSE IF(STATUS .EQ. 5) THEN
         CALL QQQSAP(CLE,DEF,VAL,N)
         CALL XIT('  CCARD ', 0 )
      ELSE IF(STATUS == 6) THEN
         PRINT *,' *** ERREUR: dans GET_COMMAND_ARGUMENT'
         CALL XIT('  CCARD ',-25)
      ELSE
         PASFINI = .FALSE.
      ENDIF

   ENDDO

   ! RETOURNER LE NOMBRE D'ARGUMENTS ASSOCIES A LA CLEF - SI DEMANDE

   IF(II .GT. 0) THEN
      II = INTERNE
   ENDIF

   RETURN

109 format(a,i4.4,a)
777 format(a,i4,a)

!  ========
CONTAINS
!  ========

!**S/P  QQQOBM - OBTENIR UN NOM DE CLEF OU UNE VALEUR D'UN ARGUMENT

   INTEGER FUNCTION QQQOBM( ARG,cltype )

      IMPLICIT NONE

      CHARACTER(LEN=8192) ARG
      integer cltype

!AUTEUR          J.CAVEEN, JANVIER 1991
!
!OBJET(QQQOBM)
!       FONCTION QUI PERMET D'ALLER CHERCHER LES ARGUMENTS D'UNE SEQUENCE
!       D'APPEL A UN PROGRAMME ET D'EN EXTRAIRE LES NOMS DE CLEFS ET LES
!       VALEURS A DONNER A CES CLEFS.
!       LA FONCTION QQQOBM RETOURNE:
!            - UNE VALEUR DE 1 SI ARG CONTIENT UN NOM DE CLEF
!            - UNE VALEUR DE 2 SI ARG CONTIENT UNE VALEUR A DONNER A UNE CLEF
!            - UNE VALEUR DE 3 SI ARG CONTIENT UN ARGUMENT POSITIONEL
!            - UNE VALEUR DE 5 SI ON DEMANDE LA SEQUENCE D'APPEL
!            - UNE VALEUR DE 0 LORSQUE TOUT EST FINI
!
!ARGUMENT:
!         ARG    SORTIE     NOM DE CLEF OU VALEUR RETOURNE
!*

      CHARACTER(LEN=8192), save :: QARGUP
      LOGICAL, save :: PASDCLE=.true., PUDCLE=.false., PREMIER=.true.
      INTEGER, save :: ARGNUM=0, QINDEX=0, INDFIN=0, NARG=0
      INTEGER I, INDEB, J, status
      character delim

      if (cltype .eq. 2) then
         delim = '='
      else
         delim = ':'
      endif
		
100   QQQOBM = 0 ; ARG = ' '

        IF(PREMIER) THEN
            NARG = COMMAND_ARGUMENT_COUNT()
            PREMIER = .FALSE.
        ENDIF

        IF(QINDEX .GE. INDFIN) THEN
        !  #  ALLER CHERCHER UN ARGUMENT
 
           ARGNUM = ARGNUM + 1
           IF (ARGNUM .GT. NARG)   RETURN
            
           QARGUP = ' '

           CALL GET_COMMAND_ARGUMENT(ARGNUM,QARGUP,INDFIN,status) 
           IF(((NARG .EQ. 1) .AND. ((INDEX(QARGUP,'-H ') .NE. 0) .OR. &
               (INDEX(QARGUP,'-h ') .NE. 0) &
               .OR. (INDEX(QARGUP,'-HELP') .NE. 0) &
               .OR. (INDEX(QARGUP,'-help') .NE. 0) ))) THEN
                QQQOBM = 5
                RETURN
            else if (status /= 0) then
                QQQOBM = 6
                RETURN
            ENDIF

            QINDEX = 1

         ENDIF
             
         ! EXTRAIRE UN NOM DE CLEF OU UNE VALEUR
                   
         IF((QINDEX .EQ. 1) .AND.  (QARGUP(1:1) .EQ. '-') .AND. (.NOT. PUDCLE)) THEN

            ! NOM DE CLEF

            PASDCLE = .FALSE.
            QQQOBM = 1
            QINDEX = 2
            IF(QARGUP(2:2) .EQ. '-') THEN
               PUDCLE = .TRUE.
               GO TO 100                   
               ! CHERCHER PROCHAIN ARG POSITIONEL
            ENDIF
         ELSE 
            ! VALEUR A DONNER
            IF(PUDCLE .OR. PASDCLE) THEN
               ! ARGUMENT POSITIONEL
               QQQOBM = 3
               ARG = QARGUP
               QINDEX = INDFIN
               RETURN
            ELSE
               QQQOBM = 2
            ENDIF
            IF((QARGUP(QINDEX:QINDEX) .EQ. delim) .OR. (QARGUP(QINDEX:QINDEX) .EQ. '=')) then
               QINDEX = QINDEX + 1
            endif
         ENDIF

         INDEB = QINDEX
         J = 1
         DO I= INDEB ,INDFIN
            IF((QARGUP(I:I) .EQ. '=') .OR. (QARGUP(I:I) .EQ. delim)) THEN
                  GOTO 500
            ELSE
                 ARG(J:J) = QARGUP(I:I)
                 QINDEX = QINDEX +1
                 J = J + 1
            ENDIF
         ENDDO

500      CONTINUE

      RETURN

   END FUNCTION QQQOBM

!**S/P  QQQOENV - OBTENIR UN NOM DE CLEF OU UNE VALEUR D'UN ARGUMENT
!                 A PARTIR D'UNE VARIABLE D'ENVIRONNEMENT CONTENANT LA
!                 SEQUENCE D'APPEL COMPLETE

   INTEGER FUNCTION QQQOENV( ARGP,CCARD_ARGS,L )

      IMPLICIT NONE

      CHARACTER(LEN=8192) ARGP, CCARD_ARGS
      INTEGER L

!AUTEUR          M. Lepine,  Fevrier 2003
!
!OBJET(QQQOENV)
!       FONCTION QUI PERMET D'ALLER CHERCHER LES ARGUMENTS D'UNE SEQUENCE
!       D'APPEL A UN PROGRAMME ET D'EN EXTRAIRE LES NOMS DE CLEFS ET LES
!       VALEURS A DONNER A CES CLEFS.
!       LA FONCTION QQQOENV RETOURNE:
!            - UNE VALEUR DE 1 SI ARG CONTIENT UN NOM DE CLEF
!            - UNE VALEUR DE 2 SI ARG CONTIENT UNE VALEUR A DONNER A UNE CLEF
!            - UNE VALEUR DE 3 SI ARG CONTIENT UN ARGUMENT POSITIONEL
!            - UNE VALEUR DE 5 SI ON DEMANDE LA SEQUENCE D'APPEL
!            - UNE VALEUR DE 0 LORSQUE TOUT EST FINI
!
!ARGUMENT:
!         ARGP    SORTIE     NOM DE CLEF OU VALEUR RETOURNE
!*
      Integer, save :: pos
      Logical, save :: debut=.true., pudcle=.false.

      Integer  i, ic, indfin
      character c, quotechar
      character(LEN=8192) arg
      integer, external :: longueur

      if (debut) then
         pos = 1
         debut = .false.
      endif

100   continue

         i = 1
         arg = ' '
         argp = ' '

         if (pos > L) then
            qqqoenv = 0        ! termine, fin de la sequence d'appel
            return
         endif

         c = ccard_args(pos:pos)
         getarg: DO                            ! obtenir la prochaine cle, valeur
            if ((pos > L) .or. (c .eq. ' ') .or. (c .eq. ':')) EXIT getarg
            if ((c .ne. '''') .and. (c .ne. '"')) then
               arg(i:i) = c
               pos = pos +1
               i = i + 1
               if (pos <= L) c = ccard_args(pos:pos)
            else
               quotechar = c
               pos = pos + 1
               c = ccard_args(pos:pos)
               quote: DO
                  if (c .eq. quotechar) EXIT quote
                  if (pos > L) then
                     print *,'CCARD, qqqoenv error: unmatched quote'
                     EXIT quote
                  endif
                  arg(i:i) = c
                  i = i + 1
                  pos = pos +1	
                  if (pos <= L) c = ccard_args(pos:pos)
               END DO quote
               pos = pos +1
               if (pos <= L) c = ccard_args(pos:pos)          
            endif
         END DO getarg

         indfin = longueur(arg)

         if ((arg(1:1) .eq. '-') .and. (.not. pudcle)) then
            qqqoenv = 1           ! nom de cle
            ic = 2                ! position de copie apres '-'
            if (arg(2:2) .eq. '-') then     ! -- passage en mode positionel
               pudcle = .true. 
               goto 100                     ! prochain argument positionel 
            endif
         else
            if (pudcle) then                ! argument positionel
               qqqoenv = 3
               ic = 1
            else
               qqqoenv = 2                  ! valeur
               ic = 1
            endif
            if (arg(1:1) .eq. '=') ic = 2   ! passer le caractere d'escape
         endif

         argp = arg(ic:indfin)
         if ((qqqoenv == 1) .and. (indfin == 2) .and. &
            ((arg(1:2) .eq. '-h') .or. (arg(1:2) .eq. '-H') .or. &
            (arg(1:5) .eq. '-help') .or. (arg(1:5) .eq. '-HELP') )) then
            qqqoenv = 5                      ! sequence d'appel demande
         endif

         pos = pos + 1          ! positionnement au debut du prochain argument

      return

   end FUNCTION QQQOENV

!**S/P - QQQTRNS - TRADUIRE ET TRANSFERER UNE VALEUR SELON LE TYPE
!
   SUBROUTINE QQQTRNS( SORTI,ENTRE,TYPE )

      IMPLICIT NONE

      CHARACTER * (*) SORTI,ENTRE
      INTEGER TYPE

!AUTEUR:   JAMES CAVEEN,  JUILLET 1991
!
!OBJET(QQQTRNS) - TRADUIRE UNE VALEUR EN MAJUSCULES/MINUSCULE SELON
!                 LE TYPE DE CLE ET TRANSFERE LE RESULTAT DANS SORTI
!
!ARGUMENTS    SORTI  SORTIE -  NOM RESULTANT DE LA TRANSFORMATION
!             ENTRE  ENTREE - NOM A TRADUIRE
!             TYPE      "   - TYPE DE TRADUCTION A APPLIQUER
!
!*
      EXTERNAL LOW2UP, UP2LOW

      IF(TYPE .EQ. 1) THEN
         CALL LOW2UP(ENTRE,SORTI)
      ELSE IF(TYPE .EQ. -1) THEN
         CALL UP2LOW(ENTRE,SORTI)
      ELSE
         SORTI = ENTRE
      ENDIF

      RETURN

   END SUBROUTINE QQQTRNS

!**S/P QQQSAP - IMPRIMER LA SEQUENCE D'APPEL AU PROGRAMME
!
   SUBROUTINE QQQSAP(CLE,DEF,VAL,N)

      IMPLICIT NONE

      INTEGER N
      CHARACTER * (*) CLE(N), DEF(N), VAL(N)

!AUTEUR:  JAMES CAVEEN, JUILLET 1991
!
!OBJET(QQQSAP) - SOUS-PROGRAMME SERVANT A IMPRIMER LA SEQUENCE D'APPEL 
!                AU PROGRAMME PRESENTEMENT EN EXECUTION.
!                QQQSAP IMPRIME SUR L'UNITE 6: 
!                   LE NOM DU PROGRAMME, LES DIFFERENTS NOM DE CLEFS
!                   ET LEUR PREMIERE ET DEUXIEME VALEUR DE DEFAUT.
!
!ARGUMENTS:
!        CLE      ENTREE  - TABLEAU CONTENANT LE NOM DES N CLEFS
!        DEF      ENTREE  - TABLEAU CONTENANT LA DEUXIEME VALEUR DE DEFAUT
!        VAL      ENTREE  - TABLEAU CONTENANT LA PREMIERE VALEUR DE DEFAUT
!                           OU LA VALEUR COURANTE DE LA CLEF
!        N        ENTREE  - NOMBRE DE CLEFS
!
      CHARACTER(LEN=8192) LENOM
      INTEGER I

      ! ON OBTIENT LE NOM DU PROGRAMME APPELANT
      CALL GET_COMMAND_ARGUMENT( 0, LENOM )

      WRITE(6,100)

      WRITE(6,*) trim( LENOM )

      DO I =1,N
         WRITE(6,200) trim( CLE(I) ),trim( VAL(I) ), trim( DEF(I) )
      ENDDO

      RETURN

100   FORMAT(' *** SEQUENCE D''APPEL ***'//)
200   FORMAT('     -',A,' [',A,':',A,']')

   END SUBROUTINE QQQSAP

END SUBROUTINE DIAG_CCARD

