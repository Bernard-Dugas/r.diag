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
      subroutine def_spectral_truncation (ncid,la,lrlmt)

      implicit none


      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'infomem.h'

      integer ncid,la,lrlmt

******
*
* AUTEUR    Guy Bergeron   mai 2003
*
*     Definir les parametres de la troncature spectrale (la et lrlmt) a 
*     partir des attributs d'une variable contenu dans le fichier netCDF.
*
* REVISIONS
*
* B. Dugas octobre 2014 :
*  - Declarations/initialisations locales des variables ktr et lmt
*  - Supprimer les variables locales nwds, lenth et jbuf
* B. Dugas, ete 2007 :
* - Allocation dynamique de memoire automatique
* - Remplacer les appels a fxlrlmt/dimgt par clrlmt/dimgt2
*
******

******CCCma

      integer*2 trunc_count
      integer   i,n, ktr,lr,lm,lmt, buf

******netCDF

      integer err
      integer varid,nlen,attnum,type,status

      character*80 name
      character*128 string

      integer, external :: clrlmt

*-----------------------------------------------------------------------

      lmt=0 ; ktr=-1
      trunc_count=lmt
      string=''

      do varid=1,nlist

         if(list(varid)%ndim .gt. 1) then

            do attnum=1,list(varid)%nattr

               status = nf_inq_attname(ncid,varid,attnum,name)
               call handle_err2(status,'def_spectral_truncation')
         
               status = nf_inq_att(ncid,varid,name,type,nlen)
               call handle_err2(status,'def_spectral_truncation')

               if (name .eq. 'trunc_count') then

                  status=nf_get_att_int2 (ncid,varid,name, trunc_count)
                  call handle_err2(status,'def_spectral_truncation')

               else if  (name .eq. 'trunc_type') then

                  status=nf_get_att_text(ncid,varid,name, string)
                  call handle_err2(status,'def_spectral_truncation')

               endif

            enddo

            if (lmt.lt.trunc_count) then

               lmt=trunc_count
               
               if (string.eq.'Triangular') ktr=2
               if (string.eq.'Rhumboidal') ktr=0
               
               lrlmt = clrlmt( lmt+1,lmt+1,ktr,.true. )
               call dimgt2(buf,la,lr,lm,ktr,lrlmt,-1)

            endif

         endif

      enddo

*-----------------------------------------------------------------------
      end
