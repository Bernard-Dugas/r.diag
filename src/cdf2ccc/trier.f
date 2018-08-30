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
      subroutine trier (NCID)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'
      include 'specmem.h'

      integer ncid

******
*
*AUTEUR Guy Bergeron       juin 2003
*
*REVISIONS
*
* B.Dugas aout 2018 :
* - Modifier les initialisations dans def_grille_lambert
* B.Dugas juil 2018 :
* - Ajouter le traitement des attributs 'inverse_flattening'
*   et 'semi_major_axis' dans la routine def_grille_lambert
* B.Dugas juil 2012 :
* - La valeur par defaut de rlonoff pour des grilles
*   tournees est de 180 degres
* - En supposant que la grille d'entree est definie de
*   (-rlonoff) a (360-rlonoff), on lui appliquera une rotation
*   de (+rlonoff). On doit alors appliquer une rotation de
*   (-rlonoff) a longpol pour lui conserver sa vrai valeur
* - dlat1=dlat2=0 => longpol=180 dans def_grille_pt,
*   donc pas de calculs pour determiner la rotation
* - Utiliser DSET_IGS plutot que ZIPIG pour
*   definir les descripteurs IG1 et IG2 des TicTacs
* - Traiter le cas 'lambert_conformal'
* - Ajouter la routine def_grille_lambert pour ce faire
* B.Dugas avr 2011 :
* - Tenir compte de minlon et maxlon lors la definition
*   des variables de reference TicTac dans def_grille_pt
* B.Dugas ete 2007 :
* - Traiter immediatement les variables connues pour les cas
*   'polar_stereographic' et 'rotated_latitude_longitude'
* - Ajouter les routines def_grille_ps et def_grille_pt
*   pour ce faire
*
******

******netCDF :

      character(256) cfield
      integer i,id,len

*-----------------------------------------------------------------------

      do id=1,nvars

      if (list(id)%var_ok) then

         list(id)%var_ok=.false.

         call up2low( list(id)%name, cfield )

*     Eliminer les variables *_bnds telles que "time_bnds" :

         if (index( cfield,'_bnds' ) /= 0) then
            if (
     .          (index( cfield,'_bnds' )+4 == 
     .           len_trim( cfield ))            
     .                      .and.
     .          (dim(list(id)%dimid(1))%name.eq.'nbnds' .or.
     .           dim(list(id)%dimid(1))%name.eq. 'bnds')
     .         ) then

               write(6,6001)trim( cfield )

            endif

*     Cas spectral, nous recuperons mean, add_offset et scale_factor :

         else if (cfield(1:4).eq.'mean') then 

            len=1
            do i=1,list(id)%ndim
               len=len*list(id)%dimid(i)
            enddo
            call get_vard(ncid,id,list(id)%type,len,mean)
            level2did=list(id)%dimid(1)
         
         else if (cfield(1:10).eq.'add_offset') then 

            len=1
            do i=1,list(id)%ndim
               len=len*list(id)%dimid(i)
            enddo
            call get_vard(ncid,id,list(id)%type,len,add_offset)

         else if (cfield(1:12).eq.'scale_factor') then 

            len=1
            do i=1,list(id)%ndim
               len=len*list(id)%dimid(i)
            enddo
            call get_vard(ncid,id,list(id)%type,len,scale_fact)


*     Nous traitons immediatement les variables caracteres connues:

         else if (list(id)%type.eq.nf_char) then

            if (cfield.eq.'polar_stereographic'        .or.
     .          cfield.eq.'rotated_latitude_longitude' .or.
     .          cfield.eq.'rotated_pole'             ) then
               call get_attribut(ncid,id,list(id)%nattr,cfield)
               if (cfield.eq.'polar_stereographic') then
                  call def_grille_ps( list(id)%nattr )
               else
                  call def_grille_pt( list(id)%nattr )
               endif
            elseif (cfield           == 'Times'       .and.
     .          dim(list(id)%dimid(1))%name == 'DateStrLen') then
               call def_times( ncid,id, dim(list(id)%dimid(1))%len,
     .                                  dim(list(id)%dimid(2))%len )
            else
               write(6,6001)trim( cfield )
            endif

*     Cas particulier de la coordonnees Lambert Conforme qui pourrait
*     ne pas etre placee dans une variable de type caractere:

         else if(cfield == 'lambert_conformal'        .or.
     .           cfield == 'lambert_conformal_conic') then

            call get_attribut(ncid,id,list(id)%nattr,cfield)
            call def_grille_lambert( list(id)%nattr )

*     Ce qui reste devrait etre une variable meteo :

         else 
            list(id)%var_ok=.true.
         endif
      endif

      enddo

      return

*-----------------------------------------------------------------------
 6001 format(/,' TRIER : neglige la variable : ',a)
*-----------------------------------------------------------------------
      end
      subroutine def_grille_ps( nattrs )

      implicit none

      integer  nattrs

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      
*     Guy Bergeron       juin 2003
*
*REVISIONS
*
* B.Dugas juillet 2008 :
* Modifier les calculs de pi, pj et dgrw selon les corrections
* apportees aux definitions des attributs false_easting, false_northing 
* et straight_vertical_longitude_from_pole dans la routine
* def_attr_polar_stereographic en juillet 2008
*
******

      real     pi      ! distance du pole selon x en nbre de dx
      real     pj      ! distance du pole selon y en nbre de dy
      real     d60     ! valeur de dx vrais a 60 degres de l'hemisphere nhem
      real     dgrw    ! angle entre l'axe des x et Greenwich(degres ouest positif)
      real     nhem    ! hemisphere nord (1), sud (2) ou global (0)
      real     latproj ! latitude_of_projection_origin (= either +90. or -90.)

      real     lamda0,stdpar,pid60,pjd60
      integer  i

*-----------------------------------------------------------------------
      nhem = -2.
      latproj = -1024.

      do i=1,nattrs
         if (attr(i)%name.eq.
     .      'straight_vertical_longitude_from_pole') then
            call attr_value( lamda0,i )
         elseif (attr(i)%name.eq.
     .      'latitude_of_projection_origin') then
            call attr_value( latproj,i )
         elseif (attr(i)%name.eq.
     .      'standard_parallel') then
            call attr_value( stdpar,i )
         elseif (attr(i)%name.eq.
     .      'false_easting') then
            call attr_value( pid60,i )
         elseif (attr(i)%name.eq.
     .      'false_northing') then
            call attr_value( pjd60,i )
         elseif (attr(i)%name.eq.
     .      'hemisphere_of_standard_parallel') then
            call attr_value( nhem,i )
         elseif (attr(i)%name.eq.
     .      'resolution_at_standard_parallel') then
            call attr_value( d60,i )
         endif            
      enddo

      if (latproj.gt.-1000.) then
         if (latproj.eq.90.) then
            nhem = 1
         else if (latproj.eq.-90.) then
            nhem = 2
         else
            print *,' latitude_of_projection_origin = ',latproj,
     .              ' non supporte'
            call xit(' def_grille_ps ',-1 )
         endif
      endif
      
      if (stdpar.ne.60.)
     .    write(6,'(/A,f10.2,A/)')
     .   ' Projection PS vraie a une latitude de ',stdpar,' degres'

      pi = 1+pid60/d60
      pj = 1+pjd60/d60

      dgrw = -(lamda0+90)
      if (dgrw < 0.  ) dgrw = dgrw+360.
      if (dgrw > 360.) dgrw = dgrw-360.

      project%len  = 8
      project%name = 'polar_stereographic'

      project%nampar(1) = 'pi'      ; project%value(1) = pi
      project%nampar(2) = 'pj'      ; project%value(2) = pj
      project%nampar(3) = 'dgrw'    ; project%value(3) = dgrw
      project%nampar(4) = 'd60'     ; project%value(4) = d60

      project%nampar(5) = 'nis'     ; project%value(5) = -1
      project%nampar(6) = 'njs'     ; project%value(6) = -1

      project%nampar(7) = 'latproj' ; project%value(7) = latproj
      project%nampar(8) = 'nhem'    ; project%value(8) = nhem

      write(6,6000) pi,pj,dgrw,d60,nhem

      return
*-----------------------------------------------------------------------

 6000 format(/22x,'pi   = ',f10.2/
     .        22x,'pj   = ',f10.2/
     .        22x,'dgrw = ',f10.2/
     .        22x,'d60  = ',f10.2/
     .        22x,'nhem = ',f10.2/)

      end
      subroutine def_grille_pt( nattrs )

      implicit none

      integer  nattrs

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'infomem.h'
      include 'varmem.h'
      include 'dimmem.h'
      include 'ztypmem.h'
      
      real dlon1 ! true longitude of the center of the computational domain
      real dlat1 ! true latitude of the center of the computational domain
      real dlon2 ! true longitude of a point on the equator of the computational domain
      real dlat2 ! true latitude of a point on the equator of the computational domain

*     CMC/RPN :

      LOGICAL           rpn_info,rpn_debug
      COMMON  /ZZVERBO/ rpn_info
      COMMON  /ZZDEBUG/          rpn_debug

      integer ZIG1,ZIG2,ZIG3,ZIG4
      real    dxla,dyla, olon,olat
      real    maxlon,minlon
      logical rotate, ok1,ok2,ok3,ok4
      integer i,j,k, ni,nj, nila,njla,xy
      real*8  theta,ro(3,3),centre(3,2)

*-----------------------------------------------------------------------
      longpol=-999.999

      ! La valeur par defaut de rlonoff pour des
      ! grilles tournees est de 180. degres

      if (rlonoff < -999.999) rlonoff = 180.

      do i=1,nattrs
         if (attr(i)%name.eq.
     .      'grid_north_pole_latitude') then
            call attr_value( gnplat,i )
         elseif (attr(i)%name.eq.
     .      'grid_north_pole_longitude') then
            call attr_value( gnplon,i )
         elseif (attr(i)%name.eq.
     .      'north_pole_grid_longitude') then
            call attr_value( longpol,i )
         endif            
      enddo

      project%len  = 5
      project%name = 'rotated_pole'

*     premier calculs des dlon1, dlat1, dlon2
*     et dlat2 en supposant que longpol=180.

      dlon1 = 180.+gnplon
      dlat2 =   0.

      if (gnplat.ge.0.0) then
         dlat1 =  90.-gnplat
         dlon2 = -90.+gnplon
      else
         dlat1 = -90.-gnplat
         dlon2 = +90.+gnplon
      endif

      ! En supposant que la grille d'entree est definie de
      ! -rlonoff a 360-rlonoff, on lui appliquera une rotation
      ! de +rlonoff. On doit alors appliquer une rotation de
      ! -rlonoff a longpol pour lui conserver sa vrai valeur

      longpol = longpol - rlonoff
      if(longpol >= 360.0)longpol=longpol-360.0
      if(longpol <    0.0)longpol=longpol+360.0

      if (dlat1 == 0. .and. dlat2 == 0.) longpol = 180.

      if (longpol.ne.180. .and. longpol.ne.-999.999) then

*        offset du pole p/r a 180

         theta=(180.-longpol)*asin(1d0)/90.

*        matrice de rotation correspondante

         ro(1,1)=+cos(theta);ro(1,2)=+sin(theta);ro(1,3)=0.
         ro(2,1)=-sin(theta);ro(2,2)=+cos(theta);ro(2,3)=0.
         ro(3,1)= 0.        ;ro(3,2)= 0.        ;ro(3,3)=1.
          
*        calculons la matrice de rotation correspondant a un
*        centre du repere tourne (dlon1,dlat1) place a 180

         CALL D_CROT( dlon1,dlat1,dlon2,dlat2 )

         if (RPN_DEBUG)
     .   write(6,6101) 'matrice rrot',((rrot(i,j),j=1,3),i=1,3)

*        re-calcul des deux premiers vecteurs de la transformation
*        ( i.e. (X) et (Y) ) avec une rotation de theta radians
*        donc [ ro x rrot ] T(-1,0,0) et [ ro x rrot ] T(0,-1,0)

         centre=0.

         do i=1,2
         do j=1,3
            do k=1,3
               centre(j,i)=centre(j,i)-ro(i,k)*rrot(k,j)
            enddo
         enddo
         enddo

         if (RPN_DEBUG)
     .   write(6,6101) 'vecteur centre',((centre(j,i),j=1,3),i=1,2)

         CALL D_CARTALL( dlon1,dlat1,centre(1,1),1 ) 
         CALL D_CARTALL( dlon2,dlat2,centre(1,2),1 ) 

      endif

      if (dlon1 .gt. 360.) dlon1=dlon1-360.
      if (dlon1 .lt.   0.) dlon1=dlon1+360.
      if (dlon2 .gt. 360.) dlon2=dlon2-360.
      if (dlon2 .lt.   0.) dlon2=dlon2+360.

      project%nampar(1) = 'dlon1' ; project%value(1) = dlon1
      project%nampar(2) = 'dlat1' ; project%value(2) = dlat1
      project%nampar(3) = 'dlon2' ; project%value(3) = dlon2
      project%nampar(4) = 'dlat2' ; project%value(4) = dlat2

      if (RPN_INFO .or. ccc_pktyp(1:2).ne.'SQ')
     .write(6,6100) dlon1,dlat1,dlon2,dlat2

      call cxgaig( 'E', ZIG1, ZIG2, ZIG3, ZIG4,
     .                  dlat1,dlon1,dlat2,dlon2 )

      ni = dim(coord(xid)%dimid(1))%len
      nj = dim(coord(yid)%dimid(1))%len

      alon(1:ni) = dcoordonne(1:ni,xid)+rlonoff
      maxlon = maxval( alon ) ; minlon = minval( alon )
!     do i=1,ni
!        if (alon(i) > 360.) alon(i) = alon(i)-360.
!        if (alon(i) <   0.) alon(i) = alon(i)+360.
      do i=2,ni
         if (alon(i) < alon(i-1)) then
            alon(i) = alon(i)+360. 
         else if (alon(i) == alon(i-1)) then
            write(6,6001) i-1,i,alon(i-1)
         endif
      enddo

*     on force la re-definition immediate de la
*     variable globale invj afin de s'assurer que
*     les valeurs de alat SOIENT croissantes

      invj = (dcoordonne(1,yid).gt.dcoordonne(nj,yid))

      if (invj) then
         do j=1,nj
            alat(nj-j+1) = dcoordonne(j,yid)
         enddo
      else
         alat(1:nj) = dcoordonne(1:nj,yid)
      endif

*     re-calcul de gnplon,gnplat et surtout de
*     rrot qui sera utilise par zipig pour
*     definir les ZIP1, ZIP2 et ZIP3 

      call d_rota( lonr,latr, alon,alat, 
     .             dlon1,dlat1,dlon2,dlat2,
     .             gnplon,gnplat, ni,nj )

      if (RPN_DEBUG .or. ccc_pktyp(1:2).ne.'SQ')
     .write(6,6102) theta*(90./asin( 1d0 )),gnplon,gnplat

      rotate = .false.
      if (gnplat.ne.90.0) rotate = .true.

*     definir les descripteur de grilles Z (>> et ^^)

      call def_lon_lat( xid, olon,dxla, nila, ok1,ok2,xy,ok3,ok4 )
      call def_lon_lat( yid, olat,dyla, njla, ok1,ok2,xy,ok3,ok4 )

      zip3 = 0 ; call dset_igs( zip1,zip2, lonr,latr,
     .                   'E',zig1,zig2,zig3,zig4,
     .                   ni,nj )

!     call zipig( ZIP1,ZIP2,ZIP3, dxla,dyla,
!    .            nila,njla, ni,nj,
!    .            rrot,rotate )

      if (ccc_pktyp(1:2).eq.'SQ') then
         call putzref( alon,alat, 'Z',
     .                 'E',ZIG1,ZIG2,ZIG3,ZIG4,
     .                 ZIP1,ZIP2,ZIP3, ni,nj )
      endif

      return
*-----------------------------------------------------------------------

 6001 format(/' Probleme avec longitudes no. ',2I5,
     .        ' = ',f10.4,'. Coordonnee non monotone.'/)

 6100 format(/21x,'dlon1 = ',f10.4/
     .        21x,'dlat1 = ',f10.4/
     .        21x,'dlon2 = ',f10.4/
     .        21x,'dlat2 = ',f10.4/)
 6101 format(/A/(10x,3d20.10))
 6102 format(/'Apres rotation de ',f10.6,' degres...'
     .       /' north_pole_grid_longitude = ',f10.4
     .       /' grid_north_pole_latitude  = ',f10.4)

      end
      subroutine def_times (ncid,varid,dim1,dim2)

      implicit none

      integer ncid,varid,dim1,dim2
 
      include 'netcdf.inc'
      include "cdf2ccc.h"
      include 'varmem.h'
      include 'infomem.h'

!     Lire/Decoder la variable caractere Times

      character (len=dim1) dates(dim2)
      integer yyyy,mm,jj,hh,mn,ss
      integer status,ii
      
      status=nf_get_var_text( ncid,varid,dates )
      call handle_err2( status,'def_times' )

      do ii=1,dim2
         read(dates(ii),1111) yyyy,mm,jj,hh,mn,ss
         dcoordonne(ii,tid) = yyyy * 1000000D0
     .                      +   mm * 10000D0
     .                      +   jj * 100D0
     .                      +   hh * 1D0
     .                      +   mn / 6D1
     .                      +   ss / 3.6D3
      enddo

 1111 format(i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2)
      
      return
      end
      subroutine def_grille_lambert( nattrs )

      implicit none

      integer  nattrs

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      
*     Bernard Dugas       juillet 2012
*
*REVISIONS
*
******

      real(8) stdpar(2) ! standard_parallel
      real(8) lonproj   ! longitude_of_central_meridian
      real(8) latproj   ! latitude_of_projection_origin
      real(8) feasting  ! false_easting
      real(8) fnorthing ! false_northing
      real(8) semajaxis ! semi_major_axis
      real(8) invflatng ! inverse_flattening

      character(256) cfield
      integer i

*-----------------------------------------------------------------------

      ! Initialiser les descripteurs

      project%name = "lambert_conformal_conic" ; project%len  = 8

      feasting  = -999999999. ; fnorthing = -999999999.
      latproj   = -999.       ; lonproj   = -999.    ; stdpar = -999.
      semajaxis = -999999999. ; invflatng = -999.

      do i=1,nattrs
         if (
     .      attr(i)%name == 'grid_mapping_name') then
            call attr_cvalue( cfield,i )
            call up2low( cfield,cfield )
            if (cfield /= "lambert_conformal_conic") return
         else if (
     .      attr(i)%name == 'standard_parallel') then
            call attr_dvalue2( stdpar(1),i,1 )
            if (attr(i)%len > 1) call attr_dvalue2( stdpar(2),i,2 )
         else if (
     .      attr(i)%name == 'latitude_of_projection_origin') then
            call attr_dvalue( latproj,i )
         else if (
     .      attr(i)%name == 'longitude_of_projection_origin' .or.
     .      attr(i)%name == 'longitude_of_central_meridian') then
            call attr_dvalue( lonproj,i )
         else if (
     .      attr(i)%name == 'false_easting') then
            call attr_dvalue( feasting,i )
         else if (
     .      attr(i)%name == 'false_northing') then
            call attr_dvalue( fnorthing,i )
         else if (
     .      attr(i)%name == 'semi_major_axis') then
            call attr_dvalue( semajaxis,i )
         else if 
     .      (attr(i)%name == 'inverse_flattening') then
            call attr_dvalue( invflatng,i )
         endif            
      enddo

      project%nampar(1) = 'stdpar1'   ; project%value(1) = stdpar(1)
      project%nampar(2) = 'stdpar2'   ; project%value(2) = stdpar(2)
      project%nampar(3) = 'latproj'   ; project%value(3) = latproj
      project%nampar(4) = 'lonproj'   ; project%value(4) = lonproj
      project%nampar(5) = 'feasting'  ; project%value(5) = feasting
      project%nampar(6) = 'fnorthing' ; project%value(6) = fnorthing
      project%nampar(7) = 'semajaxis' ; project%value(7) = semajaxis
      project%nampar(8) = 'invflatng' ; project%value(8) = invflatng


      write(6,6000) stdpar,latproj,lonproj,feasting,fnorthing,
     .              semajaxis,invflatng

      return
*-----------------------------------------------------------------------

 6000 format(/22x,'stdpar(1) = ',f15.2/
     .        22x,'stdpar(2) = ',f15.2/
     .        22x,'latproj   = ',f15.2/
     .        22x,'lonproj   = ',f15.2/
     .        22x,'feasting  = ',f15.2/
     .        22x,'fnorthing = ',f15.2/
     .        22x,'semajaxis = ',f15.2/
     .        22x,'invflatng = ',f15.2/)

      end
