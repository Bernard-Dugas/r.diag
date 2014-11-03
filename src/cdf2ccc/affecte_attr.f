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
      subroutine affecte_attr(nbr,TYPE,NAME,NLEN,
     .                           STRING,I1VALS,I2VALS,IVALS,RVALS,DVALS)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'

      integer nbr,type,nlen
      character*(*) name,string
      integer*1 i1vals(nlen)
      integer*2 i2vals(nlen)
      integer   ivals(nlen)
      real*4    rvals(nlen)
      real*8    dvals(nlen)

******
*
*AUTEUR Guy Bergeron    juin 2003
*
*     Affecte les informations relatives au type derive attr(i)
*
*REVISIONS
*
*     Bernard Dugas  fev 2009
*     Ajouter le support des donnees de type nf_byte
*
******

      integer i
*-----------------------------------------------------------------------

      nbr=nbr+1
      attr(nbr)%type=type
      attr(nbr)%name=name
      attr(nbr)%len =nlen
      
      if (type.eq.nf_char) then          ! character

         attr(nbr)%cvalue=string(1:nlen)

      else if (type.eq.nf_byte) then    ! integer*1
 
        do i=1,nlen
            attr(nbr)%i1value(i)=i1vals(i)
         enddo

      else if (type.eq.nf_short) then    ! integer*2
 
        do i=1,nlen
            attr(nbr)%i2value(i)=i2vals(i)
         enddo

      else if (type.eq.nf_int) then      ! integer
 
        do i=1,nlen
            attr(nbr)%ivalue(i)=ivals(i)
         enddo

      else if (type.eq.nf_float) then    ! real*4

         do i=1,nlen
            attr(nbr)%rvalue(i)=rvals(i)
         enddo

      else if (type.eq.nf_double) then   ! real*8

         do i=1,nlen
            attr(nbr)%dvalue(i)=dvals(i)
         enddo

      endif

*-----------------------------------------------------------------------
      end
