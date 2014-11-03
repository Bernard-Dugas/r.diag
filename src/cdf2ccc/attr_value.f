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
      subroutine attr_value2 (value,id,n)

      implicit none

      include 'cdf2ccc.h'
      include 'netcdf.inc'

      integer id,n
      real    value

******
*
*AUTEUR Guy Bergeron       avril 2004
*
*     Extraire la valeur associe a l'attribut et affecte une variable "real"
*
*REVISIONS
*
*     Bernard Dugas juillet 2012
*     - Renommer attr_value et attr_dvalue a attr_value2 et
*       attr_dvalue2 suite a l'ajout du nouvel argument n
*     - Tenir compte de attr(id)%len (qui doit etre >= n)
*     - Les routines attr_value et attr_dvalue deviennent
*       des points d'entree pour lesquels n=1
*     Bernard Dugas fevrier 2009
*     - Ajouter le support des donnees de type nf_byte
*
******

******netCDF :


******CCCma :

      LOGICAL           rpn_info,rpn_debug
      COMMON  /ZZVERBO/ rpn_info
      COMMON  /ZZDEBUG/          rpn_debug

      integer, save :: ni=1

******
*-----------------------------------------------------------------------

      ni = n

      if (attr(id)%len < n) then
         if (RPN_DEBUG)
     .      write(6,6011) id,attr(id)%type,attr(id)%len,
     .                    trim( attr(id)%name ),n
         call                                    xit('attr_value2', -1 )
      endif

      entry attr_value (value,id)

      if (attr(id)%type.eq.nf_byte) then
         
         value = real( attr(id)%i1value(ni) )
            
      else if (attr(id)%type.eq.nf_short) then
         
         value = real( attr(id)%i2value(ni) )
            
      else if (attr(id)%type.eq.nf_int) then

         value = real( attr(id)%ivalue(ni) )

      else if (attr(id)%type.eq.nf_float) then

         value = attr(id)%rvalue(ni)

      else if (attr(id)%type.eq.nf_double) then

         value = real( attr(id)%dvalue(ni) )

      else

         if (RPN_DEBUG) then
             write(6,6001) id,attr(id)%type,attr(id)%len,
     .                  trim( attr(id)%name )
             write(6,6002) nf_byte,nf_short,nf_int,nf_float,nf_double
         endif
         call                                    xit('attr_value2', -2 )

      endif

      ni = 1

*-----------------------------------------------------------------------
 6001 format('id,attr(id)%type,attr(id)%len,attr(id)%name =',3i6,A)
 6002 format(' les types connus sont...'/
     .       '  nf_byte nf_short  nf_int  nf_float nf_double '/5(i5,4x))
 6011 format('id,attr(id)%type,attr(id)%len,attr(id)%name =',3i6,A/
     .       ', requesting item no. ',i3)

      end
      subroutine attr_dvalue2(value,id,n)

      implicit none

      include 'cdf2ccc.h'
      include 'netcdf.inc'

      integer id,n
      real*8  value

*
*     Bernard Dugas      avril 2008 (lourdement inspire de attr_value)
*
*     Extraire la valeur associe a l'attribut et affecte une variable "real*8
*
*     Bernard Dugas      juillet 2012
*     Tenir compte de attr(id)%len
*
******netCDF :


******CCCma :

      LOGICAL           rpn_info,rpn_debug
      COMMON  /ZZVERBO/ rpn_info
      COMMON  /ZZDEBUG/          rpn_debug

      integer, save :: ni=1

******

*-----------------------------------------------------------------------

      ni = n

      if (attr(id)%len < n) then
         if (RPN_DEBUG)
     .      write(6,6011) id,attr(id)%type,attr(id)%len,
     .                    trim( attr(id)%name ),n
         call                                   xit('attr_dvalue2', -1 )
      endif

      entry attr_dvalue (value,id)

      if (attr(id)%type.eq.nf_byte) then
         value = dble( attr(id)%i1value(ni) )
            
      else if (attr(id)%type.eq.nf_short) then
         value = dble( attr(id)%i2value(ni) )
            
      else if (attr(id)%type.eq.nf_int) then

         value = dble( attr(id)%ivalue(ni) )

      else if (attr(id)%type.eq.nf_float) then

         value = dble( attr(id)%rvalue(ni) )

      else if (attr(id)%type.eq.nf_double) then

         value =       attr(id)%dvalue(ni)

      else

         if (RPN_DEBUG) then
             write(6,6001) id,attr(id)%type,attr(id)%len,
     .                  trim( attr(id)%name )
             write(6,6002) nf_byte,nf_short,nf_int,nf_float,nf_double
         endif
         call                                   xit('attr_dvalue2', -2 )

      endif

      ni = 1

*-----------------------------------------------------------------------
 6001 format('id,attr(id)%type,attr(id)%len,attr(id)%name =',3i6,1x,A)
 6002 format(' les types connus sont...'/
     .       '  nf_byte nf_short  nf_int  nf_float nf_double '/5(i5,4x))
 6011 format('id,attr(id)%type,attr(id)%len,attr(id)%name =',3i6,A/
     .       ', requesting item no. ',i3)

      end
      subroutine attr_cvalue (value,id)

      implicit none

      include 'cdf2ccc.h'
      include 'netcdf.inc'

      integer id
      character*(*) value

*
*     Bernard Dugas       juillet 2012
*
*     Extraire la valeur associe a l'attribut et affecte une variable "character"
*

******netCDF :
      integer ni
      character*128 string

******CCCma :

      LOGICAL           rpn_info,rpn_debug
      COMMON  /ZZVERBO/ rpn_info
      COMMON  /ZZDEBUG/          rpn_debug

******
*-----------------------------------------------------------------------

      if (attr(id)%type == nf_char) then
         
         call clean_char( attr(id)%cvalue,string,ni )
         value=string(1:ni)
  
      else

         if (RPN_DEBUG) then
             write(6,6001) id,attr(id)%type,attr(id)%len,
     .                  trim( attr(id)%name )
             write(6,6002) nf_char
         endif
         call                                        xit('attr_cval',-1)

      endif

*-----------------------------------------------------------------------
 6001 format(/'id,attr(id)%type,attr(id)%len,attr(id)%name =',3i6,A)
 6002 format(' le seul type connu est... nf_byte =',i5/)

      end
