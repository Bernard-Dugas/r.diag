!
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
!
      subroutine les_arg( CLES,DEF,def1,DEF2,NBR,VERSION )

      implicit none

      integer    nbr
      character*(*) cles(nbr),def(nbr),def1(nbr),def2(nbr),version

******
*
* AUTEUR Guy Bergeron     juilllet 2003
*
*     Decoder les cles passees en appel du programme. Recherche la cle 
*     "cles(i)";- si trouvee:  a) lire sa valeur si elle est la
*                              b) prendre le defaut (def2) si elle n'est pas la
*
*               - si non trouvee : prendre le defaut de def1
*
*     cles : les cles a cherche et decoder
*     def1 : valeur de sortie et defaut si il n'y a pas la cle
*     def2 : defaut si il y a juste la cle
*
*     NOTA: pour valeurs negatives, appeler avec -cle =-valeur
*
*     NOTA: il peut y avoir plus d'une valeur associee a un cle.
*
*           i.e.:   pgm -cle1 toto tata titi -cle2 coucou -cl3
*
*           Alors "def1" a la sortie contient toto,tata et titi.  
*     
* REVISIONS
*
* Bernard Dugas decembre 2018 :
* - Ajouter le support de '_' a la fin des cles (tel que dans CCARD).
*   La presence de ce caractere '_' minusculise la valeur du champs
* Bernard Dugas 30 janvier 2018 :
* - Remplacer les commandes F77 GETARG et IARGC par
*   GET_COMMAND_ARGUMENT et COMMAND_ARGUMENT_COUNT,
*   respectivement
* Bernard Dugas juillet 2007 :
* - La boucle principale devient un "Do while"
* - Call xit en cas d'erreur plutot que "stop"
* - Agrandir les formats I/O dans la routine definition
* Anne Frigon Oct 2003 : Adapte pour avoir cles de valeurs negatives
*
******

      integer*4     argc,iarg,jarg
      integer       i,j,m,n,init,fini,nt,nlen
      character(512) string,dummy
*-----------------------------------------------------------------------

      argc = command_argument_count()         ! nombre total d arguments d appel

      if(.false.) print*,'argc =',argc

      nt=len(string)         !longueur max permise pour string

      iarg=1

      do while (iarg <= argc) !boucle sur tous les arguments d appel

         call get_command_argument(iarg,string)

         if(string(1:1) == '-') then                ! chercher une cle

            string=string(2:nt)
            if (string == 'h')call definition(cles,def,def1,nbr,version)
            jarg=iarg+1

            n=0
            do j=1,nbr
               i=len_trim(string)
               if(cles(j) == string            .or.
     .            cles(j) == string(1:i)//'_') then ! definir la cle
                  
                  n=j ; m=len_trim(cles(n))
                  if(jarg <= argc)then              ! affecter une valeur

                     init=1
 100                 call get_command_argument(jarg,string)       ! boucler jusqu'a la cle suivante
                     if(string(1:1) == '-') then

                        if (init == 1) def1(n)=def2(n) !1ere cle sans valeur

                     else if(string(1:1) == '=')then   !possible valeur negative

                        call def_name(string,nt,' ',dummy,nlen)
                        dummy=dummy(2:nlen)       !sans le =
                        fini=init+nlen-1
                        def1(n)(init:fini)=dummy
                        init=fini+1
                        jarg=jarg+1

                     else if(string(1:1) /= ' ')then

                        call def_name2(string,nt,dummy,nlen)
                        if(init.eq.1)def1(n)=' '
                        dummy=dummy(1:nlen)
                        if(cles(n)(m:m) == '_')then
                           call up2low(dummy(1:nlen),dummy(1:nlen))
                        endif
                        fini=init+nlen
                        def1(n)(init:fini)=dummy(1:nlen)
                        init=fini+1
                        jarg=jarg+1
                        goto 100

                     endif
                  else
                     def1(n)=def2(n)     !juste la cle sans valeur
                  endif

                  iarg=jarg
                  exit

               endif
            enddo               !j

              if(n == 0) then
                 call def_name(string,nt,' ',dummy,nlen)
                 write(6,6001) string (1:nlen)
                 call xit('les_arg',-1 )
              endif
         endif

      enddo                     !i

*----------------------------------------------------------------------
 6001 format(/,'LES_ARG : "-',a ,'"'," n'est pas une cle valide !",/)
*----------------------------------------------------------------------
      end


*----------------------------------------------------------------------
*----------------------------------------------------------------------
      subroutine definition(CLES,DEF,DEF1,NBR,version)

*
* Impression pour le help
*
      implicit none

      integer    i,nbr,nblen
      character*(*) cles(nbr),def(nbr),def1(nbr),version

*----------------------------------------------------------------------

      write(6,5999)version(1:nblen(version))
      write(6,6000)
      do i=1,nbr
         write(6,6001) cles(i),def(i),def1(i)
      enddo
      write(6,6001)
      call                 xit('les_arg', 1)

*----------------------------------------------------------------------
 5999 format(/,25x,a)
 6000 format(/X,'NOMS',4x,': ',
     .       '(TYPE) DEFINITIONS',42x,' : ','DEFAUTS'/)
                 
 6001 format(x,a8,': ',a60,' : ',a)
*----------------------------------------------------------------------
      end

