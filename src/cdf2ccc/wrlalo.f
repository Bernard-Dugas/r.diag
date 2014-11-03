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
      subroutine wrlalo2( LAT_UNIT,LON_UNIT )


      implicit none

      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'
      include 'varmem.h' 
      include 'ibufmem.h'
      include 'workmem.h'

      integer lat_unit,lon_unit

******
*
*AUTEUR Guy Bergeron       avril 2004
*
*     Ecrire les valeurs des variables de coordonnees lon et lat dans des
*     fichiers de sortie.
*
*REVISIONS 
*
*     Bernard Dugas fevrier 2014 :
*     - Faire appel a COMBLINE5 (ajouter FILL_CDF_NAN=.false.)
*     Bernard Dugas juin 2013 :
*     - Faire appel a COMBLINE4 (enlever les miss_*)
*     Bernard Dugas aout 2012 :
*     - Remplacer COMBLINE2 par COMBLINE3 (deux nouveaux arguments)
*     Bernard Dugas avril 2008 :
*     - Valeur par defaut de opack = -16 (i.e. lorsque npack = 999)
*     Bernard Dugas hiver 2007 :
*     - Nouveaux arguments lat_unit,lon_unit passes a putfld2.
*       On change donc le nom de la routine a wrlalo2
*     - Allocation de memoire avec allocate/deallocate
*     - Allouer la memoire pour les grilles tournees
*     - Modifier l'appel a combline
*     G. Bergeron aout 2005 :  Reorganisation des arguments d'appels de combline
*     G. Bergeron   juin 2005: definition de invj en parametre d'appel.
*     C. Desrochers 17 mars 2005: ajouter parametres fill_cdf_def et fill_cdf
*                              a l'appel de combline ( inutilise ici)
*     Anne Frigon 14mars 2005: ajouter parametres miss_cdf_def et miss_cdf
*                              a l'appel de combline ( inutilise ici)
*
******

******netCDF :

      integer  lonid,latid
      logical  miss_all,fill_all

******CCCma :

      integer ccctime
      character*4 cccname
      integer i,j,ij,indice,ilevel,nblen
      integer nlon, nlat, nhem, opack

      real*8  scale,offset,bad,miss_cdf,fill_cdf

      data scale,offset,miss_cdf,fill_cdf /1.0,0.0,0.0,0.0/
      data ilevel,ccctime /1,0/
      data lonid,latid /0,0/
******
      real(8), dimension (:), allocatable :: dlon,dlat
*-----------------------------------------------------------------------
      allocate ( dlon(maxlen),dlat(maxlen) )

                        opack = npack
      if (opack == 999) opack = -16

*     Identifier idlon et idlat:
         
      do i=1,ncoord
         if(coord(i)%name.eq.trim(lon)) then
            lonid=i            
         else if(coord(i)%name.eq.trim(lat)) then
            latid=i 
         endif
      enddo

      if (lonid.eq.0 .or. latid.eq.0) then
         write(6,1020)
         call                                          xit('wrlalo2',-1)
      endif

*  Ouverture des fichiers CCCma

      write(6,1010)

*     Generer un tableau de lon et de lat:

      nlon=dim(coord(lonid)%dimid(1))%len
      nlat=dim(coord(latid)%dimid(1))%len

      do j=1,nlat
         do i=1,nlon
            ij=(j-1)*nlon+i
            dlon(ij)=dcoordonne(i,lonid)
            dlat(ij)=dcoordonne(j,latid)
         enddo
      enddo

           
*     Ecriture dans le fichier de sortie :

      do i=1,project%len
         if(project%nampar(i).eq.'nhem')nhem=int(project%value(i))
      enddo

*     Latitude :

      call def_cccma(coord(latid)%name,coord(latid)%mult,
     .                                     coord(latid)%add,cccname,bad)   

      indice=0
      call combline5( dlat,dval,indice,1,nlon,nlat,1,scale,
     .            offset,coord(latid)%mult,coord(latid)%add,
     .            .false.,fill_cdf,invj,0, fill_all,.false. )

      call setlab(ibuf,'GRID',ccctime,cccname,ilevel,nlon,nlat,
     .                                                       nhem,opack)

      call putfld2( lat_unit, dval, ibuf,maxpk )


*     Longitude :

      call def_cccma(coord(lonid)%name,coord(lonid)%mult,
     .                                     coord(lonid)%add,cccname,bad)   

      indice=0
      call combline5( dlon,dval,indice,1,nlon,nlat,1,scale,
     .            offset,coord(lonid)%mult,coord(lonid)%add,
     .            .false.,fill_cdf,invj,0, fill_all,.false. )

      call setlab(ibuf,'GRID',ccctime,cccname,ilevel,nlon,nlat,
     .                                                       nhem,opack)

      call putfld2( lon_unit, dval, ibuf,maxpk )

      deallocate ( dlon,dlat )
*-----------------------------------------------------------------------
 1010 format(/,/," Creation des fichiers : LONGITUDE et LATITUDE ")
 1020 format(/,/," ATTENTION :"
     ./,"     Il n'existe pas de coordonnees LONGITUDE et/ou LATITUDE",
     ./,"     dans le fichier source",/)
      end
