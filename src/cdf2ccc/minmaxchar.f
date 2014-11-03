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
      subroutine minmaxchar (IN,out,NCHAR)
      
      implicit none

      integer NCHAR
      character*(*) IN,out

******
*
*AUTEUR Guy Bergeron    juin 2003
*
*     Realise la conversion minuscule majuscule ou vice versa
*
*
* REVISIONS:
*
* G. Bergeron  juin 2005
*     Generalisation a tout caracteres alphanumerique.
*
* Anne Frigon Juillet 2004 : 
*    Ajoute lettres w et W manquantes dans data lower et upper...
*
******

      integer i

*-----------------------------------------------------------------------
      do i=1,nchar

         if(ichar(in(i:i)).ge.97.and.ichar(in(i:i)).le.122) then
            out(i:i)= char(ichar(in(i:i))-32)
         elseif(ichar(in(i:i)).ge.65.and.ichar(in(i:i)).le.90) then
            out(i:i)= char(ichar(in(i:i))+32)
         else
            out(i:i)=in(i:i)
         endif

      enddo
*-----------------------------------------------------------------------
      end


