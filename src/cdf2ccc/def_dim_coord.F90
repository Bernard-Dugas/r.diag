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
      subroutine def_dim_coord(FUNIT,IBUF)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'
      include 'varmem.h'
      include 'workmem.h'

      integer funit,ibuf(*)

!*****
!
!AUTEUR Guy Bergeron   juin 2003
!
!     Definition du monbre de dimensions (ndims), des noms dim(id)%name et 
!     de leurs longueurs dim(id)%len. Nous definissons aussi les coordonnees,
!     c'est a dire, le nom coord(id)%name, le nombre de dimension coord(id)%ndim,
!     le vecteur des ID corespondant au dimension coord(id)%dimid et les valeurs
!     de la coordonnee coordonne(i,id).
!
!REVISIONS
!
!  Guy Bergeron juillet 2004 : Hybrid height coordinate
!
!AUTRES REVISIONS
!
!  B. Dugas mai '17 : 
!  - Convertir en fichier .F90 pour traiter
!    le macro taille_entete avec s.f90
!  B.Dugas octobre '14 :
!  - Declarations/initialisations locales des variables ktr et lmt
!  B.Dugas aout '12 :
!  - Meme si nlev=1, lorsque toutes les variables partagent ce niveau
!    (i.e. invar(1)%unique_L=T), la coordonnee verticale est definie
!  B.Dugas juin '12 :
!  - Ajouter le support de 'Log Pressure Hybrid Levels' (VKIND=5002)
!  - Ne plus definir de coordonnes verticale si nlev =1
!  B.Dugas mai '12 :
!  - Definir la dimension bnds
!  - Ajouter le support de 'plev' et 'height' pour les noms
!    de la coordonnee verticale
!  B.Dugas avril '12 :
!  - Lorsque level_desc corresponds a une coordoonne qui suit le terrain,
!    coord(zid)%mult=1000. seulement si Is_CCC=T
! B. Dugas, decembre 2011:
! - Definir 'mult=1' pour 'Hybrid Height'
! B. Dugas, octobre 2008:
! - Ajouter le support de coordonnees arbitraire, niveau de sol et TOA
! B. Dugas, ete/automne/hiver 2007:
! - Ajouter le support des coordonnees de pression hybride et de hauteur
! - Ajouter le support des fichiers standards CMC/RPN.
!
!*****

      integer, parameter :: head = taille_entete

      character*8 xname,yname,zname
      character*4 data_type

      integer jbuf(head)

      integer la,lr,lm,lrlmt, ktr,lmt
      integer i,ii,ntime,nlev

      integer lablvl(maxlev)

      logical ok,Is_CCC

      data ok /.false./
      integer n
!-----------------------------------------------------------------------
      if      (trim(level_desc) == 'Pressure Levels' ) then 
         zname = 'plev'
      else if (trim(level_desc) == 'Height'            .or. &
               trim(level_desc) == '10 m'              .or. &
               trim(level_desc) == '2 m'             ) then
         zname = 'height'
      else
         zname = 'lev'
      endif

      range(1)=vlarge ; range(2)=-vlarge        ! inititaliser min et max

      Is_CCC = .false. ; if (ccc_pktyp(1:2) == 'PK') Is_CCC = .true.

!     Identifier le nbre de pas de temps, le nbre de niveau et leurs etiquettes :

      call scanfile(funit,ibuf,ntime,nlev,lablvl)

      do i=1,head
         jbuf(i)=ibuf(i)
      enddo

      if (project%name     .eq.'gaussian'             .or.  &
         (project%name(1:7).eq.'lon/lat'              .and. &
          project%name     .ne.'lon/lat regional')    .or.  &
          project%name     .eq.'polar_stereographic') then

         !Representation: global(0), hem. nord(1), hem.sud(2)
         project%nampar(project%len)='nhem'
         project%value(project%len)=float(ibuf(7))

      endif

!     Identifier le type de troncature et le nombre d'ondes (cas spectral) :

      write(data_type,'(a4)') jbuf(1) 

      if (data_type.eq.'SPEC') then

         spec=.true.
                                               ! ktr=truncation type          
         call dimgt(ival,la,lr,lm,ktr,jbuf(7)) ! la =fctn(lrlmt)=((lmt+1)**2+(lmt+1))/2
         lmt=lm-1                              ! lmt=truncation count (m)

      else                                     ! cas aux points de grilles

!     Definir les dimensions horizontales:

         ii=0

         if (project%name.eq.'rotated_latitude_longitude'  .or. &
             project%name.eq.'rotated_pole'              ) then

            xname='rlon'
            yname='rlat'

         else if (project%name(1:7).eq.'lon/lat'   .or. & 
                  project%name     .eq.'gaussian') then

            xname=trim(lon)
            yname=trim(lat)

         else

            xname='xc'
            yname='yc'

            endif

         if (.not.Is_CCC) then
!           seulement eliminer repetition de Greenwich pour grilles B,
!           les grilles G et A n'etant pas supposees en avoir besoin
            if(project%name.eq.'lon/lat global B') ii=1
         else
            if(project%name(1:14).eq.'lon/lat global' .or. &
               project%name      .eq.'gaussian'     ) then
               if(mod(jbuf(5),2).gt.0)ii=1 ! eliminer repetition de Greenwich 
            endif
         endif
         
         call affecte_dim(ndims,xdid,jbuf(5)-ii,xname)
         call affecte_dim(ndims,ydid,jbuf(6)   ,yname)
         dim(xdid)%duplic=ii            ! grille decalee (1/0) (shifted grid)

      end if

      if (spec) call affecte_dim(ndims,numdid,la*2,'num_values')

      if(level_desc /= "Surface"            .and. &
         level_desc /= "Sea Level"          .and. &
        (nlev > 1  .or. infvar(1)%unique_L)) then
         ok=.true.
         call affecte_dim(ndims,zdid,nlev,zname)
      endif

!     Definir la dimension temporelle :

      call affecte_dim(ndims,timedid,ntime,'time')
      if (time_bnds_L) &
      call affecte_dim(ndims,bndsdid,2,'bnds')

!     Definir la dimension verticale:
 
!     Definir les coordonnees horizontale :

      call affecte_coord(ncoord,xid,nf_double,1,xdid,dim(xdid)%name)
      call affecte_coord(ncoord,yid,nf_double,1,ydid,dim(ydid)%name)

      if (ok) then

!        Definir la cordonnee verticale :

         call affecte_coord( ncoord,zid,nf_double,1,zdid,trim(zname) )

!NOTA: Definir "mult" et "add" de sorte que level=value*mult+add
!      soit egale a la sortie de "lvdcode". La consistance de l'encodage
!      et du decodage est assuree dans "def_level"
!
         coord(zid)%add=0.0 ; coord(zid)%mult=1.0

         if      (trim(level_desc) == 'Sigma Levels'                .or. &
                  trim(level_desc) == 'Gal-Chen Levels'             .or. &
                  trim(level_desc) == 'Hybrid Levels'               .or. &
                  trim(level_desc) == 'Log Pressure Hybrid Levels') then
         
            if (Is_CCC) coord(zid)%mult=1000.0

         else if (trim(level_desc) == 'Pressure Levels' ) then

            coord(zid)%mult=0.01

         else if (trim(level_desc) /= 'Height'           .and. &
                  trim(level_desc) /= '10 m'             .and. &
                  trim(level_desc) /= '2 m'              .and. &
                  trim(level_desc) /= 'Arbitrary Levels' .and. &
                  trim(level_desc) /= 'Soil Layers'      .and. &
                  trim(level_desc) /= 'Top of Atmosphere'.and. &
                  trim(level_desc) /= 'Hybrid Height'  ) then

            write(6,6001)trim(level_desc)
            call                            xit('def_dim_coord' , -2)

         endif

      endif

!     Definir la coordonnee temporelle :

      call affecte_coord(ncoord,tid,nf_double,1,timedid,'time')

!     Evaluer les valeurs de la coordonnee verticale :

      if(ok) call def_level(lablvl,'decode')       


!     Evaluer les valeurs des coordonnees horizontales :

      if (project%name(1:7).eq.'lon/lat'.or. &
                                        project%name.eq.'gaussian') then

         call eval_lonlat(dim(xdid)%len,dim(ydid)%len)

      else if (project%name.eq.'polar_stereographic') then

         call eval_xcyc(dim(xdid)%len,dim(ydid)%len)

      else if (project%name.eq.'rotated_latitude_longitude'  .or. &
               project%name.eq.'rotated_pole'              ) then

         call eval_rlonlat(dim(xdid)%len,dim(ydid)%len)

      else

         write(6,6002)trim(project%name)
         call                                  xit('def_dim_coord' , -3)

      endif

!     Definir les valeurs de la coordonnee temporelle :

      do i=1,dim(timedid)%len
         dcoordonne(i,tid)=rtime(i)
      enddo

!     La dimension pour htoit

      if(level_desc.eq.'Gal-Chen Levels') &
                                  call affecte_dim(ndims,ii,1,'zt')
!-----------------------------------------------------------------------
 6001 format(' def_dim_coord : cas non definis pour level_desc = -', &
                                                                  a,'-')
 6002 format(' def_dim_coord : grille non definie pour grid_desc = -', &
                                                                  a,'-')
 6099 format(' ',A4,I10,2X,A4,I10,2I6,I9,I3)     ! format de ibuf
!-----------------------------------------------------------------------
      end subroutine def_dim_coord
