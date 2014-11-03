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
      subroutine justifie_droite(cccname)

      implicit none

      character*(*) cccname

******
*
*AUTEUR Guy Bergeron         Juin 2005
*    
*     Justifie a droite une chaine de caractere cccname.
*
******

      integer i,nt,nlen,nblen,len
      
      character*128 xname
*-----------------------------------------------------------------------
      xname=cccname

      nt=len(cccname)        ! longueur totale de cccname
      nlen=nblen(xname)      ! nombre de caracteres dans xname

      do i=1,nt
         cccname(i:i)=' '
      enddo

      i=nt-nlen+1 
      cccname(i:nt)= xname(1:nlen)
*-----------------------------------------------------------------------
      return
      end
