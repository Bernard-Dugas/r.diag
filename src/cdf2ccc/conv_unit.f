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
*-----------------------------------------------------------------------
      subroutine conv_unit(dvalue,mult,add,ni)

      implicit none
      

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'infomem.h' 

      integer i,ni
      real*8 mult,add
      real*8 dvalue(ni)

******
*
*AUTEUR Guy Bergeron             avril 2004
*
*     Effectue la conversion d'unite cccma->netcdf d'une variable.
* 
*REVISIONS
*
*     Bernard Dugas, avril 2008 : On travaille maintenant en "real*8".
*     
******

*-----------------------------------------------------------------------

      do i=1,ni
         dvalue(i)=(dvalue(i)-add)/mult
      enddo

*-----------------------------------------------------------------------
      end
