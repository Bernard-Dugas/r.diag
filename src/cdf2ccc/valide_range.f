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
      subroutine valide_range(id,valide,scale,offset)

      implicit none
      

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'infomem.h' 

      integer i,id
      integer*2 valide(max_len)

      real*8 scale,offset

******
*
*AUTEUR Guy Bergeron             avril 2004
*
*     Evalue les valeurs de valide_range a partir de unpacked_valide_range
*     et de scale (scale_factor) et offset (add_offset).
*
*REVISIONS
*
*     Bernard Dugas fevrier 2009 :
*     - Ajouter le support des donnees de type nf_byte
*
*     Bernard Dugas avril 2008 :
*     - Arguments scale,offset sont declares "real*8"
*     
******

*-----------------------------------------------------------------------

      do i=1,attr(id)%len

         if(attr(id)%type.eq.nf_double) 
     .                  valide(i)=int((attr(id)%dvalue(i)-offset)/scale)


         if(attr(id)%type.eq.nf_float) 
     .                  valide(i)=int((attr(id)%rvalue(i)-offset)/scale)


         if(attr(id)%type.eq.nf_short) 
     .                 valide(i)=int((attr(id)%i2value(i)-offset)/scale)


         if(attr(id)%type.eq.nf_byte) 
     .                 valide(i)=int((attr(id)%i1value(i)-offset)/scale)

      enddo
*-----------------------------------------------------------------------
      end
