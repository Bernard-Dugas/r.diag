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
      subroutine affecte_coord(nd,id,TYPE,NDIM,DIMIDS,NAME)

      implicit none

      include 'cdf2ccc.h'
      include 'infomem.h'

      integer nd,id,type,ndim,dimids(ndim)
      character*(*) name

******
*
*AUTEUR Guy Bergeron    juin 2003
*
*     Affecte les informations relatives au type derive coord(i)
*
******

      integer i
*-----------------------------------------------------------------------
      nd=nd+1   
      id=nd
      coord(id)%name = name
      coord(id)%type = type
      coord(id)%ndim = ndim
      do i=1,ndim
         coord(id)%dimid(i)= dimids(i)
      enddo
*-----------------------------------------------------------------------
      end
