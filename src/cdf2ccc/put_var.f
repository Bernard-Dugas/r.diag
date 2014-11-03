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
      subroutine put_var(NCID,VARID,TYPE,LEN,VALUE,SCALE,OFFSET,CONV)

      implicit none
      
      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'    

      integer ncid,varid,type,len
      real*4  value(len)
      real*8  scale,offset

      logical conv

******
*
*AUTEUR Guy Bergeron         juin 2003
*
*     Ecrire en fonction du "type" les valeurs de la variable "varid" en 
*     un bloc dans le fichier netCDF. La variable d'entre est convertie 
*     en {i2val=int((value - offset)*scal)} facultativement.
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
*     Bernard Dugas  fev 2009
*     Ajouter le support des donnees de type nf_byte
*
*     Bernard Dugas  mai 2008
*     Ajouter la routines put_vara
*
*     B. Dugas automne 2007 :
*     Syntaxe vectorielle de F90 et remplacer les INTs par des NINTs
*
******

      integer status,i

*-----------------------------------------------------------------------

      if (conv) then  ! definir i1val/i2val = int((value - offset)*scale)

         if (type == nf_byte) then
            if (scale.ne.0.0) then
               i1val(1:len) = nint( (value(1:len)-offset) /scale)
            else
               i1val(1:len) = nint( value(1)-offset )
            endif
            status=nf_put_var_int1 (ncid,varid,i1val)
         elseif (type == nf_short)then
            if (scale.ne.0.0) then
               i2val(1:len) = nint( (value(1:len)-offset) /scale)
            else
               i2val(1:len) = nint( value(1)-offset )
            endif
            status=nf_put_var_int2 (ncid,varid,i2val)
         endif

      else if (type == nf_byte)then
         do i=1,len
            i1val(i)=nint(value(i))
         enddo
         status=nf_put_var_int1(ncid,varid, i1val) 

      else if (type.eq.nf_short)then
         do i=1,len
            i2val(i)=nint(value(i))
         enddo
         status=nf_put_var_int2(ncid,varid, i2val) 

      else if (type.eq.nf_int)then
         do i=1,len
            ival(i)=nint(value(i))
         enddo
         status=nf_put_var_int(ncid,varid, ival) 

      else if (type.eq.nf_float) then

         status=nf_put_var_real(ncid,varid, value) 
         
      else if (type.eq.nf_double)then
         do i=1,len
            dval(i)=dble(value(i))
         enddo
         status=nf_put_var_double(ncid,varid, dval) 
      endif
      call handle_err2(status,'put_var')

*-----------------------------------------------------------------------
      end
      subroutine put_vara (NCID,VARID,TYPE,LEN,START,COUNT,VALUE,
     .                                                SCALE,OFFSET,CONV)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'

      integer ncid,varid,type,len
      integer start(maxdim),count(maxdim)
      real*4  value(len)
      real*8  scale,offset

      logical conv

******
*
*AUTEUR Guy Bergeron         juillet 2003
*
*     Ecrire une partie du bloc en fonction du type de la variable value.
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
******

      integer i,status

*-----------------------------------------------------------------------

      if (conv) then  ! definir i1val/i2val = int((value - offset)*scale)

         if (type == nf_byte) then

            if (scale.ne.0.0) then
               i1val(1:len) = nint( (value(1:len)-offset) /scale)
            else
               i1val(1:len) = nint( value(1)-offset )
            endif

            status=nf_put_vara_int1 (ncid,varid,start,count,i1val)

         else if (type == nf_short) then

            if (scale.ne.0.0) then
               i2val(1:len) = nint( (value(1:len)-offset) /scale)
            else
               i2val(1:len) = nint( value(1)-offset )
            endif

            status=nf_put_vara_int2 (ncid,varid,start,count,i2val)

         endif

      else if (type == nf_byte) then
         
         do i=1,len
            i1val(i)=nint(value(i))
         enddo

         status= nf_put_vara_int1(ncid,varid,start,count,i1val)
         
      else if (type.eq.nf_short) then
         
         do i=1,len
            i2val(i)=nint(value(i))
         enddo

         status= nf_put_vara_int2(ncid,varid,start,count,i2val)
         
      else if(type.eq.nf_int) then

         do i=1,len
            ival(i)=nint(value(i))
         enddo

         status= nf_put_vara_int (ncid,varid,start,count,ival)

      else if(type.eq.nf_float) then

         status= nf_put_vara_real(ncid,varid,start,count,value)

      else if(type.eq.nf_double) then

         do i=1,len
            dval(i)=dble(value(i))
         enddo

         status= nf_put_vara_double(ncid,varid,start,count,dval) 
         
      endif

      call handle_err2(status,'put_vara')

*-----------------------------------------------------------------------
         end
      subroutine put_vard2(NCID,VARID,TYPE,LEN,DVALUE,
     .                     SCALE,OFFSET,CONV,MSS,FILLVAL)

      implicit none
      
      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'    

      logical mss
      integer ncid,varid,type,len
      real*8  dvalue(len),scale,offset,fillval

      logical conv

******
*
*AUTEUR Guy Bergeron         avril 2004
*
*     Ecrir en fonction du "type" les valeurs de la variable "varid" en 
*     un bloc dans le fichier netCDF. La variable d'entree REAL*8 est 
*     convertie en {i2val=int((dvalue - offset)*scal)} facultativement.
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
* B. Dugas ete 2013 :
* - Ajouter les arguments MSS et FILLVAL (associes aux valeurs manquantes)
* - Renommer a PUT_VARD2
*
* B. Dugas automne 2007 :
* Syntaxe vectorielle de F90 et remplacer les INTs par des NINTs
*
******

      integer status,i

*-----------------------------------------------------------------------

      if (conv) then  ! definir i1val/i2val = int((value - offset)*scale)

         if (type == nf_byte) then

            if (scale.ne.0.0) then
               if (mss) then
                  do i=1,len
                     if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                        i1val(i) = int( (dvalue(i)-offset)/scale )
                     else
                        i1val(i) = nint( fillval )
                     endif
                  enddo
               else
                  i1val(1:len) = nint( (dvalue(1:len)-offset) /scale)
               endif
            else
               if (mss) then
                  do i=1,len
                     if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                        i1val(i) = nint( (dvalue(i)-offset) )
                     else
                        i1val(i) = nint( fillval )
                     endif
                  enddo
               else
                  i1val(1:len) = nint( dvalue(1:len)-offset )
               endif
            endif

            status=nf_put_var_int1 (ncid,varid,i1val)

         else if (type.eq.nf_short)then

            if (scale.ne.0.0) then
               if (mss) then
                  do i=1,len
                     if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                        i2val(i) = nint( (dvalue(i)-offset)/scale )
                     else
                        i2val(i) = nint( fillval )
                     endif
                  enddo
               else
                  i2val(1:len) = nint( (dvalue(1:len)-offset) /scale)
               endif
            else
               if (mss) then
                  do i=1,len
                     if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                        i2val(i) = nint( (dvalue(i)-offset) )
                     else
                        i2val(i) = nint( fillval )
                     endif
                  enddo
               else
                  i2val(1:len) = nint( dvalue(1:len)-offset )
               endif
            endif

            status=nf_put_var_int2 (ncid,varid,i2val)

         endif

      else if (type.eq.nf_byte)then

         if (mss) then
            do i=1,len
               if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                  i1val(i) = nint( dvalue(i) )
               else
                  i1val(i) = nint( fillval )
               endif
            enddo
         else
            i1val(1:len) = nint( dvalue(1:len) )
         endif

         status=nf_put_var_int1(ncid,varid, i1val) 

      else if (type.eq.nf_short)then

         if (mss) then
            do i=1,len
               if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                  i2val(i) = nint( dvalue(i) )
               else
                  i2val(i) = nint( fillval )
               endif
            enddo
         else
            i2val(1:len) = nint( dvalue(1:len) )
         endif

         status=nf_put_var_int2(ncid,varid, i2val) 

      else if (type.eq.nf_int)then

         if (mss) then
            do i=1,len
               if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                  ival(i) = nint( dvalue(i) )
               else
                  ival(i) = nint( fillval )
               endif
            enddo
         else
            ival(1:len) = nint( dvalue(1:len) )
         endif

         status=nf_put_var_int(ncid,varid, ival) 

      else if (type.eq.nf_float) then

         if (mss) then
            do i=1,len
               if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                  rval(i) = real( dvalue(i) )
               else
                  rval(i) = real( fillval )
               endif
            enddo
         else
            rval(1:len) = real( dvalue(1:len) )
         endif

         status=nf_put_var_real(ncid,varid, rval) 
         
      else if (type.eq.nf_double)then

         if (mss) then
            do i=1,len
               if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                  dval(i) = dvalue(i)
               else
                  dval(i) = fillval
               endif
            enddo
         else
            dval(1:len) = dvalue(1:len)
         endif

         status=nf_put_var_double(ncid,varid, dval) 

      endif

      call handle_err2(status,'put_vard2')

*-----------------------------------------------------------------------
      end


      subroutine put_varda2(NCID,VARID,TYPE,LEN,START,COUNT,DVALUE,
     .                               SCALE,OFFSET,CONV,MSS,FILLVAL)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'

      logical mss
      integer ncid,varid,type,len
      integer start(maxdim),count(maxdim)
      real*8  dvalue(len),scale,offset,fillval

      logical conv

******
*
*AUTEUR Guy Bergeron         juillet 2003
*
*     Ecrir en fonction du "type" les valeurs de la variable "varid" en 
*     partie dans le fichier netCDF. La variable d'entre est convertie 
*     en {i2val=int((value - offset)*scal)} facultativement.
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
* B. Dugas ete 2013 :
* - Ajouter les arguments MSS et FILLVAL (associes aux valeurs manquantes)
* - Renommer a PUT_VARDA2
*
* B. Dugas automne 2007 :
* Syntaxe vectorielle de F90 et remplacer les INTs par des NINTs
*
******

      integer i,status

*-----------------------------------------------------------------------

      if (conv) then  ! definir i1val/i2val = int((value - offset)*scale)

         if (type == nf_byte) then 

            if (scale.ne.0.0) then
               if (mss) then
                  do i=1,len
                     if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                        i1val(i) = int( (dvalue(i)-offset)/scale )
                     else
                        i1val(i) = nint( fillval )
                     endif
                  enddo
               else
                  i1val(1:len) = nint( (dvalue(1:len)-offset) /scale)
               endif
            else
               if (mss) then
                  do i=1,len
                     if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                        i1val(i) = nint( (dvalue(i)-offset) )
                     else
                        i1val(i) = nint( fillval )
                     endif
                  enddo
               else
                  i1val(1:len) = nint( dvalue(1:len)-offset )
               endif
            endif

            status=nf_put_vara_int1 (ncid,varid,start,count,i1val)

         else if (type.eq.nf_short) then

            if (scale.ne.0.0) then
               if (mss) then
                  do i=1,len
                     if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                        i2val(i) = nint( (dvalue(i)-offset)/scale )
                     else
                        i2val(i) = nint( fillval )
                     endif
                  enddo
               else
                  i2val(1:len) = nint( (dvalue(1:len)-offset) /scale)
               endif
            else
               if (mss) then
                  do i=1,len
                     if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                        i2val(i) = nint( (dvalue(i)-offset) )
                     else
                        i2val(i) = nint( fillval )
                     endif
                  enddo
               else
                  i2val(1:len) = nint( dvalue(1:len)-offset )
               endif
            endif

            status=nf_put_vara_int2 (ncid,varid,start,count,i2val)

         endif

      else if (type.eq.nf_byte) then
         
         if (mss) then
            do i=1,len
               if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                  i1val(i) = nint( dvalue(i) )
               else
                  i1val(i) = nint( fillval )
               endif
            enddo
         else
            i1val(1:len) = nint( dvalue(1:len) )
         endif

         status= nf_put_vara_int1(ncid,varid,start,count,i1val)
         
      else if (type.eq.nf_short) then
         
         if (mss) then
            do i=1,len
               if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                  i2val(i) = nint( dvalue(i) )
               else
                  i2val(i) = nint( fillval )
               endif
            enddo
         else
            i2val(1:len) = nint( dvalue(1:len) )
         endif

         status= nf_put_vara_int2(ncid,varid,start,count,i2val)
         
      else if(type.eq.nf_int) then

         if (mss) then
            do i=1,len
               if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                  ival(i) = nint( dvalue(i) )
               else
                  ival(i) = nint( fillval )
               endif
            enddo
         else
            ival(1:len) = nint( dvalue(1:len) )
         endif

         status= nf_put_vara_int (ncid,varid,start,count,ival)

      else if(type.eq.nf_float) then

         if (mss) then
            do i=1,len
               if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                  rval(i) = real( dvalue(i) )
               else
                  rval(i) = real( fillval )
               endif
            enddo
         else
            rval(1:len) = real( dvalue(1:len) )
         endif

         status= nf_put_vara_real(ncid,varid,start,count,rval)

      else if(type.eq.nf_double) then

         if (mss) then
            do i=1,len
               if (abs( dvalue(i)-fill_ccc ) >= fill_toler) then
                  dval(i) = dvalue(i)
               else
                  dval(i) = fillval
               endif
            enddo
         else
            dval(1:len) = dvalue(1:len)
         endif

         status= nf_put_vara_double(ncid,varid,start,count,dval) 
         
      endif
      call handle_err2(status,'put_varda2')

*-----------------------------------------------------------------------
         end





