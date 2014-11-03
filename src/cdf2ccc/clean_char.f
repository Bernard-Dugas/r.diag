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
      subroutine clean_char(string,chaine,nlen)

      implicit none

      integer nlen
      character*(*) string,chaine

******
*
*AUTEUR Guy Bergeron             avril 2004
*
*     Elimine le caractere special '\n' d'une chaine de caracteres 
*
*Revisions
*
*  Bernard Dugas mars 2008 : Enlever les char(0) "nul" en fin de chaine
*
******     

      integer   niin,niout,newline,ib
      character(len=1024) dummy

      data newline /10/
*-----------------------------------------------------------------------

      
      niin=len_trim(string) ; niout=len(chaine)
      call get_string(string,niin,achar(newline),dummy,nlen)
      nlen=min(nlen,niout) ; chaine=dummy(1:nlen)
      do ib=1,nlen
         if (chaine(ib:ib).eq.achar(0)) then
            nlen=ib-1
            exit
         endif
      enddo
      chaine = chaine(1:nlen)

*-----------------------------------------------------------------------
      end
            
