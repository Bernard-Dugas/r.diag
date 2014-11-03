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
      subroutine def_name (string,ni,delim,name,nlen)

******
*
* AUTEUR Guy Bergeron   juin 2003
*
*     Extrait la chaine de caractere precedant le delimiteur "delim"
*     et determine sa longueur. Nous eliminons tout les blancs contenuent
*     dans la chaine sauf si celle-ci commence par des appostrophes ou
*     des guillemets
*            
* REVISIONS
*
*     B. Dugas, ete 2007 :
*     - Ajouter le point d'entree def_name2 afin de
*       tenir compte des guillemets simples et doubles
*
******

      integer       nlen,long,iga
      character*(*) name,string,delim
      character     ga

*-----------------------------------------------------------------------
      name=''
      nlen=0

      do i=1,ni
         if (string(i:i) .ne. ' ') then
            if (string(i:i) .ne. delim) then
               nlen=nlen+1
               name(nlen:nlen)=string(i:i)
            else
               return
            endif
         endif
      enddo

      return
*-----------------------------------------------------------------------

      entry def_name2 (string,ni,name,nlen)

                              ga=" "
      if (string(1:1).eq.'"') ga='"'
      if (string(1:1).eq."'") ga="'"

      long=len_trim( string )

      if (ga.ne.' ') then

         iga=index( string(2:long),ga )

         if (iga.eq.long-1) then
            name=string(2:long-1)
            nlen=long-2
         else
            call xit('def_name',-1 )
         endif

      else

         name=string(1:long)
         nlen=long

      endif

      return

*-----------------------------------------------------------------------
      end



