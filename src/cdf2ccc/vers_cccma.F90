#     if !defined (taille_entete)
#         define   taille_entete 32
#     endif
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
      subroutine vers_cccma2 (ncid,FUNIT,LAT_UNIT,LON_UNIT)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'workmem.h'
      include 'dimmem.h'
      include 'infomem.h' 
      include 'varmem.h'
      include 'specmem.h'
      include 'ibufmem.h'
      include 'ztypmem.h'

      integer ncid,funit,lat_unit,lon_unit


!*****
!
!AUTEUR Guy Bergeron         juillet  2003
!
!
!     Traduction : netCDF -> CCCma
!
!
!REVISIONS 
!
!  Bernard Dugas mai 2017 : 
!  - Convertir en fichier .F90 pour traiter
!    le macro taille_entete avec s.f90
!  Bernard Dugas octobre 2014 :
!  - Corriger la sequence d'appel a def_spectral_truncation
!  Bernard Dugas juillet 2013 :
!  - Deplacer l'appel a MISPAR vers CDF2CCC (main)
!  Bernard Dugas mai 2012 :
!  - Introduit le support de 'time_bnds'
!  Bernard Dugas novembre 2009 :
!  - Correction pour max1d > maxlen, donc
!    petites grilles et beaucoup d'echantillons
!  Bernard Dugas octobre 2008 :
!  Initialiser jbuf == 0 (requis sous IRIX64)
!  Bernard Dugas hiver 2007 :
!  - Enlever l'include de 'machtyp.h'.
!  - Nouveaux arguments lat_unit,lon_unit passe a wrlalo.
!    Donc, le nom de la routine devient vers_cccma2
!  - Allocation de memoire avec allocate/deallocate
!  - Allouer la memoire pour les grilles tournees
!  - On n'utilise plus main_memory.h
!  - Initialiser ibuf a zero au depart
!  Guy Bergeron ??? 200? : initialisation de maxpk
!
!*****

      integer, parameter :: head = taille_entete

!*****CCCma

      integer nwds,length,jbuf(head),la,lrlmt,nhem,opack

!*****netCDF

      integer ii,i,id,status,ngatts
      integer ndim,nvar,ngatt,unlimd

!-----------------------------------------------------------------------
                        opack = npack
      if (opack.eq.999) opack =-64

      jbuf = 0

!     Ouverture du fichier netCDF:

      status=nf_open(netcdf_file,nf_nowrite,ncid)
      call handle_err2(status,'vers_cccma2')

      status=nf_inq(ncid,ndims,nvars,ngatt,unlimdimid)
      call handle_err2(status,'vers_cccma2')

!     Allocation de memoire (dim,coord,var,list):

      maxdim=ndims +4
      maxvar=nvars +3
!
      allocate( dim  (maxdim) , coord(maxdim) , &
                var  (maxvar) , list (maxvar) )

!     Les parametres dimensionnelles :

      call init_dim

      call get_dim(ncid)

      max1d=1
      do i=1,ndims
         if(max1d.lt.dim(i)%len)max1d=dim(i)%len
      enddo      

!     Etablir la liste des variables du fichier netCDF :

      call define_list(ncid,nvars)


!     Identifier les variables coordonnees :

      call get_coord2(ncid)

!     Longueur maximale

      if (spec) then

         maxlen=dim(numdid)%len*dim(coord(zid)%dimid(1))%len

         call def_spectral_truncation (ncid,la,lrlmt)

         call setlab(jbuf,'SPEC',0,'TOTO',1,la,1,lrlmt,opack)

      else

         call test_dim ()
            
         if (coord(xid)%dimid(1).eq.-1 .or. &
             coord(yid)%dimid(1).eq.-1) then
            write(6,6001)
            call xit('vers_cccma2',-1 )
         endif

         maxlen=(dim(coord(xid)%dimid(1))%len+1) &
                *dim(coord(yid)%dimid(1))%len    &
                *dim(coord(zid)%dimid(1))%len

         do i=1,project%len
            if(project%nampar(i).eq.'nhem')nhem=int(project%value(i))
         enddo

         call setlab(jbuf,'GRID',0,'TOTO',1, &
            dim(coord(xid)%dimid(1))%len+1,dim(coord(yid)%dimid(1))%len, &
                                                             nhem,opack)

      endif

      maxlen=max(max1d,maxlen)
      maxlev=dim(coord(zid)%dimid(1))%len
      maxtime=dim(coord(tid)%dimid(1))%len

      call lblchk(length,nwds,opack,jbuf)
      maxpk=length-head                    !longueur (maxpk) du buffer (ibuf)

      allocate ( variable   (  maxlen,maxvar ) , &
                 dcoordonne (  max1d ,maxdim ) , &
                 time_bnds  (  2     ,max1d  ) )

      allocate ( ibuf       (  maxlen+head ), &
                 i1val      (  maxlen ) ,     &
                 i2val      (  maxlen ) ,     &
                 ival       (  maxlen ) ,     &
                 dval       (2*maxlen ) ,     &
                 rtime      (  maxlen ) ,     &
                 rval       (  maxlen ) )

      allocate(  alon       (  maxlen ) , &
                 alat       (  maxlen ) , &
                 lonr       (  maxlen ) , &
                 latr       (  maxlen ) )

      allocate ( add_offset (  maxlev*maxtime ) , &
                 scale_fact (  maxlev*maxtime ) , &
                 mean       (2*maxlev*maxtime ) )

      ibuf = 0
      ibuf(1:head) = jbuf(1:head)

!     Lire les valeurs de variables coordonnees :

      call get_coordonne(ncid)

!     Faire le tri dans les variables :

      call trier (ncid)

      if (spec) then

         call rdspec2(ncid,funit)
      else

         if (lalo) call wrlalo2( lat_unit,lon_unit )     ! Ecrire les longitudes et latitudes en sortie
         call rdlatlon2(ncid,funit)
      end if

!     Relacher la memoire

      deallocate( add_offset,scale_fact,mean,            &
                  ibuf,i1val,i2val,ival,dval,rtime,rval, &
                  list,var,coord,dim )

      deallocate( alon,alat,lonr,latr )

      return
!-----------------------------------------------------------------------
 6001 format(/' Definition manquante des coordonnees X et/ou Y'/)
      end
