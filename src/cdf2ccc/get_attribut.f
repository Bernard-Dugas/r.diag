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
      subroutine get_attribut(NCID,VARID,NATTS,VARNAME)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'

      integer ncid,varid,natts
      character*(*) varname

******
*
* AUTEUR Guy Bergeron             mai 2003
*
*     Extraire les attributs de la variables VARID dans le fichier netCDF.
*
* REVISIONS
*
*     Bernard Dugas, fevrier 2009 :
*     - Ajouter le support des donnees de type nf_byte
*     Bernard Dugas fevrier 2008 :
*     - Verifier la validite des contenus de la variable level_desc
*       et de l'attribut du meme nom avant de definir la variable
*     - La variable project%name peut maintenant etre definie
*       par les attributs 'grid_mapping' et 'grid_desc' en
*       plus de l'attribut 'grid_mapping_name'
*
******     

      integer*1 i1vals
      integer*2 ndum,trunc_count,i2vals(max_len)

      integer   i,n,ni,nlen,attnum,type,outype,status
      integer   ivals(max_len)
      real*4    rdum,rvals(max_len)
      real*8    dvals(max_len)
      logical   ok1,ok2

      character*1 check
      character*80 name
      character*128 string,string1,string2
*****
*-----------------------------------------------------------------------

      nlen=len_trim(varname)
      write(6,6000) varname(1:nlen)

      n=0
      do attnum=1,natts

         status = nf_inq_attname(ncid,varid,attnum,name)
         call handle_err2(status,'get_attribut')

         status = nf_inq_att(ncid,varid,name,type,nlen)
         call handle_err2(status,'get_attribut')

         if (type.eq.nf_char) then

            status=nf_get_att_text(ncid,varid,name, string)

         else if (type.eq.nf_byte) then

            status=nf_get_att_int1 (ncid,varid,name, i1vals)
            
         else if (type.eq.nf_short) then

            status=nf_get_att_int2 (ncid,varid,name, i2vals)
            
         else if (type.eq.nf_float) then

            status=nf_get_att_real (ncid,varid,name, rvals)

         else if (type.eq.nf_double) then

            status=nf_get_att_double (ncid,varid,name, dvals)

         endif

         call handle_err2(status,'get_attribut')
         call affecte_attr(n,type,name,nlen,
     .                       string,i1vals,i2vals,ivals,rvals,dvals)

         if (attr(n)%name.eq.'level_desc') then

            call clean_char(attr(n)%cvalue,string,ni)     !Eliminer "\n"

            if (level_desc == ' ') then
               level_desc=string
            else if (string .ne. level_desc) then
               ok1=.false. ; ok2=.false.
               do i=1,nlvl
                  if(level_desc.eq.possible%level(i))ok1 = .true.
                  if(string    .eq.possible%level(i))ok2 = .true.
               enddo
               if (.not.ok1 .or. ok2) then
                  write(6,6100) trim( level_desc ),trim( string )
                  level_desc=string
               endif
            endif
         endif

         if (attr(n)%name.eq.'grid_mapping_name' .or.
     .       attr(n)%name.eq.'grid_mapping'      .or.
     .       attr(n)%name.eq.'grid_desc'       ) then

            call clean_char(attr(n)%cvalue,string,ni)     !Eliminer "\n"
            project%name=string

         endif

         if (attr(n)%type.eq.nf_char) then
        write(6,*)attr(n)%name(1:25),' = ',attr(n)%cvalue(1:nlen)
      else if (attr(n)%type.eq.nf_byte) then
        write(6,*)attr(n)%name(1:25),' = ',(attr(n)%i1value(i),i=1,nlen)
      else if (attr(n)%type.eq.nf_short) then
        write(6,*)attr(n)%name(1:25),' = ',(attr(n)%i2value(i),i=1,nlen)
      else if (attr(n)%type.eq.nf_int) then
        write(6,*)attr(n)%name(1:25),' = ',(attr(n)%ivalue(i),i=1,nlen)
      else if (attr(n)%type.eq.nf_float) then
        write(6,*)attr(n)%name(1:25),' = ',(attr(n)%rvalue(i),i=1,nlen)
      else if (attr(n)%type.eq.nf_double) then
        write(6,6200)' '//attr(n)%name(1:25)//' = ',
     .                   (attr(n)%dvalue(i),i=1,nlen)
      endif

      enddo
      nattr=n
*-----------------------------------------------------------------------
 6000 format(/' ATTRIBUTS (',a,')')
 6100 format(/' Remplacement du type de niveau ',A,' par ',A)
 6200 format(A,(3d17.10))
*-----------------------------------------------------------------------
      end
