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
      subroutine get_var(NCID,ID,TYPE,LEN,value)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'    

      integer ncid,id,type,len
      real*4  value(len)


******
*
*AUTEUR Guy Bergeron   juin 2003
*
*     Extrait le bloc au complet en fonction du type de la variable value.
*
*     -------------------------------------------------------
*     |       |                      |             |        |
*     | XTYPE | netCDF/CDL Data Type | FORTRAN     |  Bits  |
*     |       |                      |             |        |
*     |-------|----------------------|-------------|--------|
*     |       |                      |             |        |
*     |   1   |  byte                | integer*1   |   8    |
*     |       |                      |             |        |
*     |   2   |  char                | character   |   8    |
*     |       |                      |             |        |
*     |   3   |  short               | integer*2   |   16   |
*     |       |                      |             |        |
*     |   4   |  int                 | integer     |   32   |
*     |       |                      |             |        |
*     |   5   |  float               | real*4      |   32   |
*     |       |                      |             |        |
*     |   6   |  double              | real*8      |   64   |
*     |       |                      |             |        |
*     -------------------------------------------------------
*
*REVISIONS
*
*     Bernard Dugas  mars 2012
*     Toujours transfer LEN items (BugFix nf_float)
*
*     Bernard Dugas  fev 2009
*     Ajouter le support des donnees de type nf_byte
*
*     Bernard Dugas  mai 2008
*     Ajouter les routines get_vara, get_vard et get_varda
*
******

      integer i,status
******DEBUG
      integer j,ni
*-----------------------------------------------------------------------
      if (type.eq.nf_byte)then
         status=nf_get_var_int1(ncid,id, i1val)
         call handle_err2(status,'get_var')
         value(1:len)=real( i1val(1:len) )

      else if (type.eq.nf_short)then
         status=nf_get_var_int2(ncid,id, i2val)
         call handle_err2(status,'get_var')
         value(1:len)=real( i2val(1:len) )

      else if (type.eq.nf_int) then
         status=nf_get_var_int(ncid,id, ival)
         call handle_err2(status,'get_var')
         value(1:len)=real( ival(1:len) )

      else if (type.eq.nf_float)then
         status=nf_get_var_real(ncid,id, rval)
         call handle_err2(status,'get_var')
         call infnan2( rval,len, .false. )
         value(1:len)=rval(1:len)

      else if (type.eq.nf_double) then
         status=nf_get_var_double(ncid,id, dval)
         call handle_err2(status,'get_var')
         call infnan2( dval,len, .true. )
         value(1:len)=real( dval(1:len) )

      endif 
*-----------------------------------------------------------------------
      end
      subroutine get_vara(NCID,VARID,TYPE,START,COUNT,LEN,value)

      implicit none
      
      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'    

      integer ncid,varid,type,len
      integer start(maxdim),count(maxdim)
      real*4  value(len)

*     Bernard Dugas  mars 2012
*     Toujours transfer LEN items (BugFix nf_float)
*
*     Extrait une partie du bloc en fonction du type de la variable value.
*
*     -------------------------------------------------------
*     |       |                      |             |        |
*     | XTYPE | netCDF/CDL Data Type | FORTRAN     |  Bits  |
*     |       |                      |             |        |
*     |-------|----------------------|-------------|--------|
*     |       |                      |             |        |
*     |   1   |  byte                | integer*1   |   8    |
*     |       |                      |             |        |
*     |   2   |  char                | character   |   8    |
*     |       |                      |             |        |
*     |   3   |  short               | integer*2   |   16   |
*     |       |                      |             |        |
*     |   4   |  int                 | integer     |   32   |
*     |       |                      |             |        |
*     |   5   |  float               | real*4      |   32   |
*     |       |                      |             |        |
*     |   6   |  double              | real*8      |   64   |
*     |       |                      |             |        |
*     -------------------------------------------------------

      integer i,status,err

*-----------------------------------------------------------------------
      
      if (type.eq.nf_byte) then     ! integer*1
         status=nf_get_vara_int1(ncid,varid,start,count,i1val)
         call handle_err2(status,'get_vara')
         value(1:len)=real( i1val(1:len) )

      else if (type.eq.nf_short) then ! integer*2
         status=nf_get_vara_int2(ncid,varid,start,count,i2val)
         call handle_err2(status,'get_vara')
         value(1:len)=real( i2val(1:len) )

      else if (type.eq.nf_int) then  ! integer*4
         status=nf_get_vara_int(ncid,varid,start,count,ival)
         call handle_err2(status,'get_vara') 
         value(1:len)=real( ival(1:len) )
        
      else if (type.eq.nf_float) then  ! real*4
         status=nf_get_vara_real(ncid,varid,start,count, rval)
         call handle_err2(status,'get_vara')
         call infnan2( rval,len, .false. )
         value(1:len)=rval(1:len)
         
      else if (type.eq.nf_double) then  ! real*8
         status=nf_get_vara_double(ncid,varid,start,count, dval)
         call handle_err2(status,'get_vara')
         call infnan2( dval,len, .true. )
         value(1:len)=real( dval(1:len) )

      end if
*
*-----------------------------------------------------------------------
      end
      subroutine get_vard(NCID,ID,TYPE,LEN,dvalue)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'    

      integer ncid,id,type,len
      real*8  dvalue(len)

*     Bernard Dugas  mars 2012
*     Toujours transfer LEN items (BugFix nf_double)
*
*     Extrait le bloc au complet en fonction du type de la variable et
*     affecte un REAL*8 dvalue par la valeur.
*
*     -------------------------------------------------------
*     |       |                      |             |        |
*     | XTYPE | netCDF/CDL Data Type | FORTRAN     |  Bits  |
*     |       |                      |             |        |
*     |-------|----------------------|-------------|--------|
*     |       |                      |             |        |
*     |   1   |  byte                | integer*1   |   8    |
*     |       |                      |             |        |
*     |   2   |  char                | character   |   8    |
*     |       |                      |             |        |
*     |   3   |  short               | integer*2   |   16   |
*     |       |                      |             |        |
*     |   4   |  int                 | integer     |   32   |
*     |       |                      |             |        |
*     |   5   |  float               | real*4      |   32   |
*     |       |                      |             |        |
*     |   6   |  double              | real*8      |   64   |
*     |       |                      |             |        |
*     -------------------------------------------------------

      integer i,status
*-----------------------------------------------------------------------
      if (type.eq.nf_byte) then
         status=nf_get_var_int1(ncid,id, i1val)
         call handle_err2(status,'get_vard')
         dvalue(1:len)=dble( i1val(1:len) )

      else if (type.eq.nf_short)then
         status=nf_get_var_int2(ncid,id, i2val)
         call handle_err2(status,'get_vard')
         dvalue(1:len)=dble( i2val(1:len) )

      else if (type.eq.nf_int) then
         status=nf_get_var_int(ncid,id, ival)
         call handle_err2(status,'get_vard')
         dvalue(1:len)=dble( ival(1:len) )

      else if (type.eq.nf_float)then
         status=nf_get_var_real(ncid,id, rval)
         call handle_err2(status,'get_vard')
         call infnan2( rval,len, .false. )
         dvalue(1:len)=dble( rval(1:len) )

      else if (type.eq.nf_double) then
         status=nf_get_var_double(ncid,id, dval)
         call handle_err2(status,'get_vard')
         call infnan2( dval,len, .true. )
         dvalue(1:len)=dval(1:len)
      endif 
*-----------------------------------------------------------------------
      end
      subroutine get_varda(NCID,VARID,TYPE,START,COUNT,LEN,dvalue)

      implicit none
      
      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'    

      integer ncid,varid,type,len
      integer start(maxdim),count(maxdim)
      real*8  dvalue(len)

******
*
*     Extrait une partie du bloc en fonction du type de la variable dvalue.
*
*     -------------------------------------------------------
*     |       |                      |             |        |
*     | XTYPE | netCDF/CDL Data Type | FORTRAN     |  Bits  |
*     |       |                      |             |        |
*     |-------|----------------------|-------------|--------|
*     |       |                      |             |        |
*     |   1   |  byte                | integer*1   |   8    |
*     |       |                      |             |        |
*     |   2   |  char                | character   |   8    |
*     |       |                      |             |        |
*     |   3   |  short               | integer*2   |   16   |
*     |       |                      |             |        |
*     |   4   |  int                 | integer     |   32   |
*     |       |                      |             |        |
*     |   5   |  float               | real*4      |   32   |
*     |       |                      |             |        |
*     |   6   |  double              | real*8      |   64   |
*     |       |                      |             |        |
*     -------------------------------------------------------
*
*REVISIONS
*
*     Bernard Dugas  mars 2012
*     Toujours transfer LEN items (BugFix nf_double)
*
******

      integer i,status,err

*-----------------------------------------------------------------------
      
      if (type.eq.nf_byte) then       ! integer*1
         status=nf_get_vara_int1(ncid,varid,start,count,i1val)
         call handle_err2(status,'get_varda')
         dvalue(1:len)=dble( i1val(1:len) )

      else if (type.eq.nf_short) then ! integer*2
         status=nf_get_vara_int2(ncid,varid,start,count,i2val)
         call handle_err2(status,'get_varda')
         dvalue(1:len)=dble( i2val(1:len) )

      else if (type.eq.nf_int) then  ! integer*4
         status=nf_get_vara_int(ncid,varid,start,count,ival)
         call handle_err2(status,'get_varda') 
         dvalue(1:len)=dble( ival(1:len) )
        
      else if (type.eq.nf_float) then  ! real*4
         status=nf_get_vara_real(ncid,varid,start,count, rval)
         call handle_err2(status,'get_varda')
         call infnan2( rval,len, .false. )
         dvalue(1:len)=dble( rval(1:len) )
         
      else if (type.eq.nf_double) then  ! real*8
         status=nf_get_vara_double(ncid,varid,start,count, dval)
         call handle_err2(status,'get_varda')
         call infnan2( dval,len, .true. )
         dvalue(1:len)=dval(1:len)

      end if
*
*-----------------------------------------------------------------------
      end
