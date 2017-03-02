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
      subroutine rdlatlon2 (NCID,FUNIT)

      use          diag_convert_ip123

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'
      include 'varmem.h' 
      include 'ibufmem.h'
      include 'workmem.h'
      include 'ztypmem.h'

      integer ncid,funit

******
*
*     Guy Bergeron       juin 2003
*
*     Traitement des variables definies sur une grille
*
*REVISIONS
*
* B.Dugas janvier '23 :
* - Utiliser PUTSAMPLZ pour sauver le nombre d'echantillons
*   svsm dans IBUF apres avoir place HIVAL et LOVAL dans la
*   section haute de IBUF lorsque leur taille est specifie 
*   avec l'argument -dt. Sinon, on sauve tim2-tim1 dans IP2.
* B.Dugas avril '16 :
* - Ajouter le support des grilles de type Y
* B.Dugas mars '15 :
* - Faire appel a COMBLINE6 (ajouter opack3 a la liste des arguments)
* - Verifier l'ordre des dimensions em X, Y et Z.
* B.Dugas fevrier '14 :
* - Utiliser la fonction IDNAN pour savoir si les attributs
*   _FillValue ou missing_value ont pour valeur NaN
* - Si tel est le cas, definir fill_ccc a sa valeur par
*   defaut, si cette variable n'a pas deja ete definie
* - De plus, desactiver les arrets dans InfNaN (appele dans GET_VARD)
* - Enfin, appeller COMBLINE5 qui tient compte de cette situation
* B.Dugas juin '13 :
* - Modifier le traitement des attributs _FillValue et
*   missing_value, donnant une totale preseance au premier
* - Ajouter la variable fill_toler (tolerance d'erreur) 
* - Ne plus consider les variables du type miss_ccc*
* - Faire appel a COMBLINE4
* B.Dugas octobre '12 :
* - Tenir compte des grilles Lat/Lon Regionnales
*   avec des latitudes inversees
* - Appeller def_level si 0 < zdid < ndims et nlev == 1
* B.Dugas septembre '12 :
* - Support explicte des grilles inconnues => 'DATA'
*   definies via project%name = "unknown"
* B.Dugas aout '12 :
* - Ne pas tenir compte des time_bnds lorsque l'attribut
*   cell_method n'exite pas ou vaut "time: point"
* - Faire appel a COMBLINE3 qui indique maintenant si le
*   champs courant est completement manquant (MISS_ALL) ou
*   bien complement remplace (FILL_ALL). Dans les deux cas,
*   ne plus sauver le champs en question puisqu'il est vide
* B.Dugas juin et juillet '12 :
* - La valeur par defaut de -rlonoff de -1000 est traduite
*   a 0.0 pour des grilles non tournees qu'on code en type Z.
*   Si dcoordonne(1,xid) < 0, on effectue une rotation pour
*   forcer les longitudes dans un interval de [0,360]
* - Ajouter la routine decode_single_level qui est utilisee
*   pour decoder la valeur de l'attribut 'actual_level' si present
* - Traiter le cas 'lambert_conformal_conic' => GRTYP='!', IG1=gribcode
* - Definir des descripteurs de grilles CMC/RPN par defaut lorsque
*   les coordonnees horizontales n'ont pas ete identifiees et/ou
*   definies (donc si xdid=0 ou ydid=0 => GRTYP=B)
* - De plus, si (x/y/z)did < 1 et coord((x/y/z)id) est definie (i.e.
*   dimid(1) > 0), on forcera (x/y/z)did = coord((x/y/z)id)%dimid(1)
* - Definir IP3 lorsque time_bnds=T (CMC/RPN)
* - Tenir compte des time_bnds pour une variable seulement si
*   l'attribut cell_method existe pour celle-ci et vaut soit
*   'time: mean' ou 'time: sum'. Ignorer tous les autres cas
* B.Dugas mai 2012 :
* - Tenir compte de la variable time_bnds et des arguments
*   '-dt' et '-dateo' lorsque direction = 'rpncmc' pour definir
*   les descripteurs temporels (DATEO,DEET,NPAS) de ces fichiers
* B.Dugas avr '12 :
* - Definir GRTYP='G' lorsque project%name='gaussian'
* - Lorsque dxcons=dycons=T et xglb/=0, tenter de determiner si
*   olon,dlon,olat,dlat correspondent a une grille de type 'L'.
*   Sinon, coder le tout en grille de type 'Z' regionale mais
*   non tournee. Utliser la nouvelle fonction logique CHECKRES.
* - Lorsque dxcons=T, mais dycons=F, on tente encore une fois
*   de coder le tout en grille de type 'Z' regionale non tournee
* - Valeur par defaut de l'etiquette definie lorsque la section
*   globale existe mais ne contient rien de pertinent la-dessus
* - Definir xdid, ydid et zdid lorsque ca n'a pas ete fait
*   auparavent (i.e. pas de variables correspondantes). On
*   choisit alors d'ecrire une grille de type 'B', si les
*   variables xid et yid sont definies
* - Mieux definir dimZid, dimYid et dimXid
* B.Dugas mai '09 :
*   Coordonnees NetCDF "[xyzt]coord" specifiees en arguments ?
* B.Dugas fev 2009 :
*   Ajouter le support des donnees de type nf_byte
* B.Dugas oct '08 :
* - Ne plus re-definir xdid,ydid,zdid (c'est deja fait dans GET_COORD2)
* - Ajouter le support de donnees non-geographiques et/ou intemporelles
* B.Dugas oct '07 a avr '08 :
* - Ne plus utiliser HPALLOC/HPDEALLC
* - Si npack = 999, utiliser var(nn)%type pour definir opack2
* - Tenir compte des IP1 dans IBUF(4) pour les fichiers CMC/RPN
* - Examiner les latitudes et longitudes pour tenter d'identifier
*   le type de grilles en presence, quitte a corriger project%name
*   et ce en utilisant la nouvelle routine def_lon_lat. Definir
*   les descripteurs de grilles CMC/RPN en consequence
* - Ne pas compacter les sorties CMC/RPN lorsqu'on est en presence
*   de donnees manquantes (miss) ou de remplissages (fill)
* - Si miss_val_cdf est vrai, miss_ccc_def le sera aussi, et dans
*   ce cas, si miss_ccc est indefini, il sera egal a miss_cdf
* - Si fill_val_cdf est vrai, fill_ccc_def le sera aussi, et dans
*   ce cas, si fill_ccc est indefini, il sera egal a fill_cdf
* B.Dugas septembre 2007 :
*   Ajouter le support des fichiers RPN/CMC
* A.Frigon aout 2006 :  
*   Eliminer nhem impose a 1 pour generaliser grille
*   polar_stereographic sur hemispheres NORD(1) SUD(2)
*   si direction=cccma
*   de pair avec cle cle_nhem ajoutee dans lire_arg.f
*   (defaut cle_nhem=99).
*   Cette nouvelle variable cle_nhem a ete introduite dans cdf2ccc.h
*   et ajoutee dans common card.
*   Tout ceci parce que trier.f neglige d'identifier
*   la variable polar_stereographic car definie de type char.
* G.Bergeron   aout    2005 :  Reorganisation des arguments d'appels de combline
* G.Bergeron   juin    2005 : parametre d'appel pour invj
* C.Desrochers mars    2005 : Ajout de l'option de remplissage "fill_ccc"
* G.Bergeron   avril   2004 : Declaration de coordonne en REAL*8
* G.Bergeron   juillet 2004 : Affectation de timeid fctn de "time"
*
******

******netCDF :

      logical  do_time_bnds
      character(len=512), dimension(:),allocatable :: att_name
      character(len=512)  cfield,timename

      integer  dimXid,dimYid,dimZid,timeid,nl
      integer  id,len,lena,natts,na,ntime,status

      integer  start(maxdim),count(maxdim)

      real(8)  scale,offset,miss_cdf,fill_cdf

******CCCma :

      character(len=4) cccname,type
      integer i,j,ii,iii,nn,k,kk,nlev,indice,nhem
      integer itime,ccctime, dim1,dim2

      logical fill_ccc_def0
      real    bad

      integer ilevel(maxlev)

******RPN/CMC :

      character(len=1)   GRTYP,pGRTYP,NULS
      character(len=128) ETIKET,string,level_value
      character(len=4),  save :: Gtyp='GRID',Ztyp='ZONL',Dtyp='DATA'

      integer, parameter :: TURBO16=134,IEEE=5,STD=1

      integer ::   ZIG1, ZIG2, ZIG3, ZIG4, idum,
     .              IG1,  IG2,  IG3,  IG4, datyp,
     .             pIG1, pIG2, pIG3, pIG4, ni,nj

      real    ::   olat  , olon  , dlat  , dlon  ,
     .             pi    , pj    , d60   , dgrw  ,
     .             dlat1 , dlon1 , dlat2 , dlon2

      logical ::   ok_lon, ok_lat, grille_Y, grille_Z, fill_cdf_nan

******general :

      integer,save :: svsm=-1,rkind=-1
      real,save    :: hival=-1.,loval=-1.

      real    ::   hold
      integer ::   xid2d, yid2d
      integer ::   npas,deet,datei,dateo,ip2,ip3
      real(8) ::   tdelta,tim1,tim2,zero8=0.0_8
      real(8) ::   dum81,dum82 ! Dummies for MISPAR

      real, allocatable, dimension(:,:) :: zlon,zlat

      CHARACTER(len=4), external :: GETYP
      logical,     external :: checkres,idnan,monotone
      external     puthic,puthigh,putsamplz,Diag_CONVIP_plus

      LOGICAL, dimension(nvars) :: fill_message

      LOGICAL           rpn_info,rpn_debug,noabort,holdl
      COMMON  /ZZVERBO/ rpn_info
      COMMON  /ZZDEBUG/          rpn_debug
      COMMON  /ZZABORT/                    noabort

******

      logical ::   fill_all, transpose_xy
      integer ::   fill_count=0, indice1,indice2
      integer ::   xglb,yglb,xyxy, iig2, opack1,opack2,opack3
      logical ::   ok,xbgrd,ybgrd,dxcons,dycons,ygauss,lxyxy
      logical ::   xincr,yincr,miss_val_cdf,fill_val_cdf
      character(1) cloche

*-----------------------------------------------------------------------
      if (ccc_pktyp(1:2) /= 'SQ' .or. ! => Fichier destination en format CCCma
     .    project%name == "unknown" .or. ! => Grille de type non supporte/inconnu
     .   (xdid == 0 .and. coord(xid)%nattr == -1)  .or. ! => Coordonnees horizontales non
     .   (ydid == 0 .and. coord(yid)%nattr == -1)) then ! identifiees dans le fichier source
         dxcons = .true. ; dycons = .false. ; ygauss = .true. ; ii = 0
         if (ccc_pktyp(1:2) /= 'SQ') goto 200
         goto 100
      else if (project%name == "lambert_conformal_conic") then
         ii = 0
         goto 100
      endif

      xglb = -1 ; yglb = -1 ; ii = 0 ; iig2 = 0
      xincr = .false. ; yincr  = .false.
      xbgrd = .false. ; ybgrd  = .false.
      dxcons = .true. ; dycons = .false. ; ygauss = .false.

      cloche=char(7)

      if (tid > 0 .and. ladate == -1) call decodate( ncid,zero8,ladate )
                  
      fill_message = .true.

      nn=0
      nhem=-1

      do i=1,project%len
         if(project%nampar(i).eq.'nhem')nhem=nint(project%value(i))
      enddo

      if (project%name.eq.'polar_stereographic' .and. nhem.lt.0)then
         if (cle_nhem.eq.99)then
            write(6,6001)
            call                                   xit('rdlatlon2', -1 )
         else
            nhem=cle_nhem
         endif
      endif

      do id=1,nvars

*        retirer certaines information des variables 1D :
         if (list(id)%ndim.eq.1 ) then
            if (list(id)%name.eq.'x'   .or.
     .          list(id)%name.eq.'x_2' .or.
     .          list(id)%name.eq.'lon' .or.
     .          list(id)%name.eq.'LON' .or.
     .          list(id)%name.eq.xcoord.or.
     .          list(id)%name.eq.'longitude') then
               call def_lon_lat( xid, olon,dlon,xyxy,
     .                           dxcons,xincr,xglb,xbgrd,lxyxy )
            elseif (list(id)%name.eq.'y'   .or.
     .              list(id)%name.eq.'y_2' .or.
     .              list(id)%name.eq.'lat' .or.
     .              list(id)%name.eq.'LAT' .or.
     .              list(id)%name.eq.ycoord.or.
     .              list(id)%name.eq.'latitude') then
               call def_lon_lat( yid, olat,dlat,xyxy,
     .                           dycons,yincr,yglb,ybgrd,ygauss)
            endif
            if (list(id)%var_ok .and.
     .         (list(id)%name /= 'ap' .and. list(id)%name /= 'b' ))
     .          list(id)%var_ok=.false.
         endif

      enddo

*     dcoordonne(:,xid) et dcoordonne(:,yid) ont peut-etre ete
*     definies dans get_coordonne meme si ces coordonnees ne
*     sont pas dans la liste de variables

      if (coord(xid)%nattr == -2)
     .   call def_lon_lat( xid, olon,dlon,xyxy,
     .                     dxcons,xincr,xglb,xbgrd,lxyxy )

      if (coord(yid)%nattr == -2)
     .   call def_lon_lat( yid, olat,dlat,xyxy,
     .                     dycons,yincr,yglb,ybgrd,ygauss)

      if (project%name     .eq.'gaussian' .or.
     .    project%name(1:7).eq.'lon/lat') then

*        verifier le type de grilles: gaussiennes, lat-lon A, B ou L ???

         grille_Y = .false. ; grille_Z = .false.


         if (dxcons) then
            if (ygauss .and. xglb.eq.0) then

               if (.not.yincr)
     .         iig2               = 1

               nhem               = yglb
               project%name       = 'gaussian'
               GRTYP              = 'G'

            elseif (dycons ) then

               write(6,6060)cloche

               if (xglb.eq.0 .and.
     .             yglb.ge.0) then

                  if (.not.yincr)
     .            iig2            = 1
                  nhem            = yglb
                  if (ybgrd) then
                     project%name = 'lon/lat global B'
                  else
                     project%name = 'lon/lat global'
                  endif

               else

                  if (cdf2_mode.eq.'cdf2rpn' .and. .NOT.yincr) then
                     invj  = .NOT.invj
                     yincr = .true.
                  endif

                  ! olon,dlon,olat,dlat peuvent-ils etre codes en grille L ?
                  ok_lon = checkres( olon,dlon )
                  ok_lat = checkres( olat,dlat )

                  if (ok_lon .and. ok_lat) then ! Oui

                     project%name    = 'lon/lat regional'

                  else ! Non ils ne le peuvent pas. On les code en grille Z

                     grille_Z = .true.

                  endif

               endif

            else

*              on a trouve une grille non-globale a resolution
*              uniforme en X et a resolution non-uniforme en Y
*              qui pourra etre decrite par une grille etiree
*              de type Z, non tournee

               grille_Z = .true.
               write(6,6061)cloche

            endif

         else if (.not.dxcons) then

*           commencer par chercher des variables lon
*           ET lat 2D afin de definir une grille Y ?

            xid2d = -1 ; yid2d = -1

            do id=1,nvars

               if (list(id)%ndim == 2 ) then
                  if (list(id)%name == 'lon' .or.
     .                list(id)%name == 'longitude') then
                     xid2d = id
                     cycle
                  elseif 
     .               (list(id)%name == 'lat' .or.
     .                list(id)%name == 'latitude') then
                     yid2d = id
                     cycle
                  endif
               endif

            enddo

            if (xid2d > 0 .and. yid2d > 0) then

               grille_Y = .true.
               write(6,6065)cloche

            else

*              cette grille devrait egalement etre decrite
*              comme une grille etiree de type Z, non tournee

               grille_Z = .true.
               write(6,6064)cloche

            end if

         endif

         if (grille_Y .or. grille_Z) then

            project%name = 'rotated_pole' 

            ni = dim(coord(xid)%dimid(1))%len
            nj = dim(coord(yid)%dimid(1))%len

            allocate( zlon(ni,nj),zlat(ni,nj) )

            if (rlonoff < -999.999) rlonoff = 0.0

            if (dcoordonne(1,xid) < -0.001_8)
     +          rlonoff = -dcoordonne(1,xid) 

            alon(1:ni) = dcoordonne(1:ni,xid)+rlonoff

            invj = (dcoordonne(1,yid) > dcoordonne(nj,yid))

            if (invj) then
               do j=1,nj
                  alat(nj-j+1) = dcoordonne(j,yid)
               enddo
            else
               alat(1:nj) = dcoordonne(1:nj,yid)
            endif
                     
            do j=1,nj ; zlon(:,j) = alon(:) ; enddo
            do i=1,ni ; zlat(i,:) = alat(:) ; enddo
          
            ZIP3 = 0 ; call dset_igs( ZIP1,ZIP2,zlon,zlat,
     +                                ZTYP,ZIG1,ZIG2,ZIG3,ZIG4,
     +                                ni,nj )

            if (grille_Z) then

               ! On ecrit une grille de type Z non-tournee

               dlat1 = 0.0 ; dlon1 = 180.-rlonoff ! ces parametres correspondent a
               dlat2 = 0.0 ; dlon2 = 270.-rlonoff ! une grille non-tournee usuelle

               if (dlon1 > 360.) dlon1 = dlon1-360.
               if (dlon1 <   0.) dlon1 = dlon1+360.
               if (dlon2 > 360.) dlon2 = dlon2-360.
               if (dlon2 <   0.) dlon2 = dlon2+360.

               call cxgaig( 'E', ZIG1, ZIG2, ZIG3, ZIG4,
     +                           dlat1,dlon1,dlat2,dlon2 )

               call putzref( alon,alat,'Z',
     +                       'E',ZIG1,ZIG2,ZIG3,ZIG4,
     +                           ZIP1,ZIP2,ZIP3,ni,nj )

            else if (grille_Y) then

               ! grille Y (collection de points)

               invj = .false.

               call get_var( ncid,xid2d,list(xid2d)%type,
     .                       dim(list(xid2d)%dimid(1))%len*
     .                       dim(list(xid2d)%dimid(2))%len,
     .                       zlon )
               call get_var( ncid,yid2d,list(yid2d)%type,
     .                       dim(list(yid2d)%dimid(1))%len*
     .                       dim(list(yid2d)%dimid(2))%len,
     .                       zlat )

               ok_lon = monotone( zlon, ni,nj, 1, xincr )
               ok_lat = monotone( zlat, ni,nj, 2, yincr )

               if (ok_lon .and. ok_lat) then ! (dlon1,dlat1) = min values -1.0
                  if (xincr) then
                      dlon1 = minval( zlon(1 , :) )
                  else
                      dlon1 = minval( zlon(ni, :) )
                  end if
                  if (yincr) then
                      dlat1 = minval( zlat(: , 1) )
                  else
                      dlat1 = minval( zlat(: ,nj) )
                  end if
                  if (xglb == -1) dlon1 = dlon1-1.0 ! This seems to be a non-global
                  if (yglb == -1) dlat1 = dlat1-1.0 ! grid that can be shifted
                  dlat2 = .10 ; dlon2 = .10
               else
                  dlat1 = 0.0 ; dlon1 = 0.0 ! selon la documentation
                  dlat2 = 1.0 ; dlon2 = 1.0 ! en ligne des grilles Y
               endif

               if (dlon1 > 360.) dlon1 = dlon1-360.
               if (dlon1 <   0.) dlon1 = dlon1+360.

               call cxgaig( 'L', ZIG1, ZIG2, ZIG3, ZIG4,
     +                           dlat1,dlon1,dlat2,dlon2 )

               zlon = merge( zlon-360.,zlon,(zlon > 360.) )
               zlon = merge( zlon+360.,zlon,(zlon <   0.) )
                  
               ZIP3 = 0 ; call dset_igs( ZIP1,ZIP2,zlon,zlat,
     +                                   ZTYP,ZIG1,ZIG2,ZIG3,ZIG4,
     +                                   ni,nj )

               call putzref( zlon,zlat,'Y',
     +                       'L',ZIG1,ZIG2,ZIG3,ZIG4,
     +                           ZIP1,ZIP2,ZIP3,ni,nj )

            endif

            deallocate( zlon,zlat )

         end if

         if (ccc_pktyp(1:2).ne.'SQ' .and. .NOT.yincr)
     .      write(6,6063) .NOT.invj,cloche

      endif

      ! tenter mettre le titre ou un autre bout
      ! de texte pertinent dans l'etiquette RPN/CMC

  100 status=nf_inq_varnatts( ncid,nf_global,natts )

      if (status == nf_noerr .and. natts > 0) then

         allocate( att_name(natts) )

         att_name(:) = ' '

         do na=1,natts
            status=nf_inq_attname( ncid,nf_global,na,att_name(na) )
            if (status .ne. nf_noerr) exit
         enddo

         ! chercher l'attribut global 'title'

         do na=1,natts
            call clean_char( att_name(na),cfield,len )
            call up2low( cfield, cfield )
            if (index(cfield,'title') > 0) exit
         enddo

         if (na > natts) then

            ! chercher l'attribut global 'experiment'

            do na=1,natts
               call clean_char( att_name(na),cfield,len )
               call up2low( cfield, cfield )
               if (index(cfield,'experiment') > 0) exit
            enddo

            if (na > natts) then

               ! chercher l'attribut global 'project'

               do na=1,natts
                  call clean_char( att_name(na),cfield,len )
                  call up2low( cfield, cfield )
                  if (index(cfield,'project') > 0) exit
               enddo

               if (na > natts) then

                  ! chercher l'attribut global 'institution'

                  do na=1,natts
                     call clean_char( att_name(na),cfield,len )
                     call up2low( cfield, cfield )
                     if (index(cfield,'institution') > 0) exit
                  enddo

               endif

            endif

         endif

         if (na <= natts) then
            status=nf_get_att_text( ncid,nf_global,
     .                              att_name(na), cfield )
            call clean_char( cfield,etiket,len )
         else
            etiket='Netcdf2CCC'
         endif

      else

         etiket='Netcdf2CCC'

      endif

*     definir les descripteur de grilles RPN/CMC

      if (project%name.eq.'polar_stereographic') then

         if (ccc_pktyp(1:2).eq.'SQ')
     .   Gtyp = 'SUBA'

         do i=1,project%len
            if(project%nampar(i).eq.'pi'  )pi  =project%value(i)
            if(project%nampar(i).eq.'pj'  )pj  =project%value(i)
            if(project%nampar(i).eq.'d60' )d60 =project%value(i)
            if(project%nampar(i).eq.'dgrw')dgrw=project%value(i)
         enddo

         if (nhem.eq.1) GRTYP = 'N'
         if (nhem.eq.2) GRTYP = 'S'

         call cxgaig( GRTYP, IG1, IG2, IG3, IG4,
     .                       pi,  pj,  d60, dgrw )

      elseif (project%name.eq.'rotated_latitude_longitude'  .or.
     .        project%name.eq.'rotated_pole'              ) then

         if (ccc_pktyp(1:2).eq.'SQ') Gtyp  = 'SUBA'

         if (.not.grille_Y) then
            GRTYP = 'Z'
         else
            GRTYP = 'Y'
         endif

         IG1 = ZIP1
         IG2 = ZIP2
         IG3 = ZIP3
         IG4 =  0

      elseif (project%name(1:7).eq.'lon/lat') then

         if (project%name.eq.'lon/lat regional') then

            if (ccc_pktyp(1:2).eq.'SQ')
     .      Gtyp  = 'SUBA'
            GRTYP = 'L'

            CALL CXGAIG( GRTYP, IG1, IG2, IG3, IG4 ,
     .                          olat,olon,dlat,dlon)

            if (ccc_pktyp(1:2).ne.'SQ')
     .      write(6,6101) olat,olon,dlat,dlon

         else if (ccc_pktyp(1:2).eq.'SQ') then

            if (project%name.eq.'gaussian'        )
     .                   GRTYP = 'G'
            if (project%name.eq.'lon/lat global'  ) then
                         GRTYP = 'A'
            else
                         GRTYP = 'B'
            endif

            if (xbgrd         .and.
     .         (GRTYP.eq.'G'  .or.
     .          GRTYP.eq.'A')) ii=-1
            if (.not.xbgrd    .and.
     .          GRTYP.eq.'B')  ii=+1

            IG1 = nhem ; IG2 = iig2 ; IG3 = 0 ; IG4 = 0

         endif

      else if (project%name == "lambert_conformal_conic") then

         Gtyp = 'SUBA' ; GRTYP = '!' ; nhem = 0
         IG1 = gribcode ; IG2 = 0 ; IG3 = 0 ; IG4 = 0

      else if (project%name == "unknown") then

         Gtyp = Dtyp ; GRTYP='X'
         IG1=0 ; IG2=0 ; IG3=0 ; IG4=0

      endif

*     Lorsque xid,yid ou zid ne sont pas dans la liste
*     des nvars variables, xdid,ydid ou zdid ne sont
*     peut-etre pas encore definis

      if (xdid < 1 .or. ydid < 1) then

*        Projection par defaut lorsque les valeurs
*        des coordonnees en x et y sont manquantes

         GRTYP='B' ; Gtyp='GRID' ; nhem = 0
         IG1=0 ; IG2=0 ; IG3=0 ; IG4=0

      endif

  200 if(xdid < 1)then
         if(coord(xid)%dimid(1) > 0)then
            xdid = coord(xid)%dimid(1)
         else
            write(6,6003) 'X'
            call                                   xit('rdlatlon2', -3 )
         endif
      endif
      if(ydid < 1)then
         if(coord(yid)%dimid(1) > 0)then
            ydid = coord(yid)%dimid(1)
         else
            write(6,6003) 'Y'
            call                                   xit('rdlatlon2', -3 )
         endif
      endif
      if(zdid < 1)then
         if(coord(zid)%dimid(1) > 0)then
            zdid = coord(zid)%dimid(1)
         else
            write(6,6003) 'Z'
            call                                   xit('rdlatlon2', -3 )
         endif
      endif

      pGRTYP=GRTYP ; pIG1=IG1 ; pIG2=IG2 ; pIG3=IG3 ; pIG4=IG4
 
      fill_ccc_def0 = fill_ccc_def

      time_bnds_L = .false.
      if (tid > 0) then
         do id=1,nvars          ! Chercher/Lire time_bnds
            if (list(id)%name == 'time_bnds') then
               len = 2*dim(tid)%len
               call get_vard(ncid,id,list(id)%type,len,time_bnds)
               time_bnds_L = .true.
            endif
         enddo
      endif

      write(6,6000)                                                   

      TRAITER_VARIABLES: do id=1,nvars

         if (list(id)%var_ok ) then

            list(id)%var_ok=.false.

            nn=nn+1

*           Definir le type derive var(nn)

            call affecte_var(nn,list(id)%name,list(id)%type,
     .                     list(id)%ndim,list(id)%dimid,list(id)%nattr)

*           Initialisation :

            dimXid=-1
            dimYid=-1
            dimZid=-1
            timeid=-1

*           Le traitement depend des dimensions

            lena=1 ; timename=' ' ; nlev=dim(zdid)%len

            do i=1,var(nn)%ndim
               start(i)=1             
               count(i)=dim(var(nn)%dimid(i))%len
               call clean_char( dim(var(nn)%dimid(i))%name,cfield,nl )
               call up2low( cfield, cfield )
               if(var(nn)%dimid(i) == unlimdimid  .or.
     .            cfield(1:nl)     == tcoord      .or.
     .            cfield(1:nl)     == 'time'      .or.
     .            cfield(1:nl)     == 't'       ) then
                  timename=cfield(1:nl)
                  count(i)=1
                  timeid=i
                  ntime=dim(var(nn)%dimid(i))%len
               endif
               if (dimZid < 0 .and.
     .             var(nn)%dimid(i) == coord(zid)%dimid(1)) dimZid=i
               if (dimYid < 0 .and.
     .             var(nn)%dimid(i) == coord(yid)%dimid(1)) dimYid=i
               if (dimXid < 0 .and.
     .             var(nn)%dimid(i) == coord(xid)%dimid(1)) dimXid=i
               if (dim(var(nn)%dimid(i))%name /= 'bnds')
     .         lena=lena*dim(var(nn)%dimid(i))%len
            enddo

            if (timeid == -1) ntime = 1

*           Lire les attributs de la variable :
            
            call get_attribut(ncid,id,var(nn)%nattr,var(nn)%name)
            call def_cccma(var(nn)%name,var(nn)%mult,var(nn)%add,
     .                                                      cccname,bad)   

            if (ccc_pktyp(1:2) == 'SQ') then
               GRTYP=pGRTYP ; IG1=pIG1 ; IG2=pIG2 ; IG3=pIG3 ; IG4=pIG4
               call fill_high( ibuf, GRTYP,IG1,IG2,IG3,IG4,etiket )
            endif

*           definir scale, offset, miss_cdf et do_time_bnds :

            scale=1.0
            offset=0.0
            miss_cdf=0.0
            fill_cdf=0.0
            level_value = ' '
            miss_val_cdf=.false.
            fill_val_cdf=.false.
            fill_cdf_nan=.false.
            do_time_bnds=.false.
            transpose_xy=.false.

            do i=1,var(nn)%nattr

               if(attr(i)%name.eq.'scale_factor')
     .                                         call attr_dvalue(scale,i)
               if(attr(i)%name.eq.'add_offset' )
     .                                        call attr_dvalue(offset,i)

               if(attr(i)%name.eq.'missing_value') then
                  call attr_dvalue(miss_cdf,i)
                  miss_val_cdf=.true.
               endif

               if(attr(i)%name.eq.'_FillValue') then
                  call attr_dvalue(fill_cdf,i)
                  if (idnan( fill_cdf,.true. )) fill_cdf_nan = .true.
                  fill_val_cdf=.true.
               endif

               if(attr(i)%name == 'cell_methods') then
                  call attr_cvalue(string,i)
                  call up2low( string,string )
                  if(string == trim( timename )//': mean'  .or.
     .               string == trim( timename )//': sum'   .or.
     .               string /= trim( timename )//': point')
     .               do_time_bnds = time_bnds_L
               endif

               if(attr(i)%name   == 'actual_level' .and.
     .            ccc_pktyp(1:2) == 'SQ'          ) then
                  call attr_cvalue(string,i)
                  call up2low( string,level_value )
               endif

            enddo

            if (fill_val_cdf .neqv. miss_val_cdf) then
               if (miss_val_cdf) then
                  if (idnan( miss_cdf,.true. )) fill_cdf_nan = .true.
                  fill_val_cdf = .true.
                  fill_cdf = miss_cdf
               endif
            endif
            
            if (fill_ccc_def0 .neqv. fill_ccc_def) then
               fill_ccc_def = fill_ccc_def0
               if (.not.fill_ccc_def) call UNSET_MISPAR( )

            endif

*           Si les attr 'missing_value' et '_FillValue'
*           n'existent pas dans le fichier cdf, on force
*           la cle de remplacement a .false.

            if ( fill_ccc_def .AND. (.NOT. fill_val_cdf) )  then
               fill_ccc_def=.false.
               if (rpn_info) write(6,6050)cloche
            endif

*           Si l'attribut '_FillValue' est present dans
*           le fichier cdf et qu'une cle correspondante
*           pour les fichiers CCC/RPN n'a pas ete definie,
*           on force cette cle a .true. et
*           on s'assure que fill_ccc=fill_cdf

            if ( (.NOT. fill_ccc_def) .AND. fill_val_cdf )  then
               if (fill_cdf_nan)                            then
                  call MISPAR( fill_ccc, dum81,dum82 ) ! Valeur par defaut
                  fill_toler = abs( fill_cdf ) * 0.001
               else
                  fill_ccc   =      fill_cdf
               endif
               fill_ccc_def  = .true.
            endif

            if (fill_ccc_def) call SET_MISPAR( fill_ccc, 0.001_8 )

*           Lire les valeurs de la variable :             

            if (dimZid > 0 .or.
     .         (nlev == 1  .and. level_value == ' '   .and.
     .          zdid >  0  .and. zdid  <=  ndims    )) then
               call def_level(ilevel,'encode') ! definir les ibuf(4)
            else
               nlev=1 ; ilevel(1)=1 ! Valeur par defaut
               if (level_value /= ' ')
     .            call decode_single_level( ilevel(1),level_value )
            endif

            len=1
            if (dimXid > 0) len=len*dim(xdid)%len
            if (dimYid > 0) len=len*dim(ydid)%len
            if (dimZid > 0) len=len*dim(zdid)%len

*           Dimensions horizontales (dim1,dim2) ?

            if (dimXid > 0 .and. dimYid > 0) then

               type=Gtyp

               iii=ii
               dim1=dim(xdid)%len
               dim2=dim(ydid)%len

               if (dimZid > 0 .and.
     .            (dimZid < dimYid .or. dimZid < dimXid)) then                
                  write(6,6004) dimXid,dimYid,dimZid
                  call                             xit('rdlatlon2', -4 )
               else if (dimXid == dimYid+1) then
                  transpose_xy=.true.
               else if (dimXid+1 /= dimYid) then
                  write(6,6004) dimXid,dimYid,dimZid
                  call                             xit('rdlatlon2', -4 )
               endif

            elseif (dimXid < 0 .and. dimYid < 0) then

               if (non_geographique) then

                  write(6,6080) 'convertie (DATA)...'

                  type=Dtyp ; GRTYP='X'
                  IG1=0 ; IG2=0 ; IG3=0 ; IG4=0

                  iii=0
                  dim1=lena/ntime
                  dim2=     ntime

                  len=dim1*dim2

                  timeid=-1 ; ntime=1
                  nlev=1 ; call convpr( ilevel(1),1.0,3,+1 )

                  if (index( list(id)%name,'_bnds' ) /= 0) then
                     if (
     .                   (index( list(id)%name,'_bnds' )+4 == 
     .                    len_trim( list(id)%name ))            
     .                             .and.
     .                   (dim(list(id)%dimid(1))%name.eq.'nbnds' .or.
     .                    dim(list(id)%dimid(1))%name.eq. 'bnds')
     .                  ) then
                       nlev=2 ; call convpr( ilevel(2),2.0,3,+1 )
                    endif
                  endif

               else

                  write(6,6080) 'ignoree...'

                  cycle TRAITER_VARIABLES

               endif

            else

*              Une seule dimension horizontale...

               type=Ztyp

               if (dimXid < 0) then
                  dim1=dim(ydid)%len
                  iii=0
               endif

               dim2=1

            endif

            if (ntime*len .ne. lena) then
               write(6,6002) ntime, len, lena
               if (dimXid > 0) then
                  write(6,6012) '1',dim(xdid)%len
               else
                  write(6,6022) '1'
               endif
               if (dimXid > 0) then
                  write(6,6012) '2',dim(ydid)%len
               else
                  write(6,6022) '2'
               endif
               if (dimXid > 0) then
                  write(6,6012) '3',dim(zdid)%len
               else
                  write(6,6022) '3'
               endif
               write(6,6020)
               call                                xit('rdlatlon2', -2 )
            endif

            if (var(nn)%type == nf_byte)   opack3 = -8
            if (var(nn)%type == nf_short)  opack3 = -16
            if (var(nn)%type == nf_float)  opack3 = -32
            if (var(nn)%type == nf_double) opack3 = -64

            if (npack == 999 .or. ntime == 1) then
               opack1 = opack3
            else
               opack1 = npack
            endif

            write(6,6020)
            do itime=1,ntime                           ! boucle temporelle

               if (tid > 0) then
                  tim1 = dcoordonne(itime,tid)
                  if (direction == 'rpncmc') then
                     if (do_time_bnds) then
                        tim1 = time_bnds(1,itime)
                        tim2 = time_bnds(2,itime)
                        if (ladate /= -1 .and. nint( dt ) >= 1) then
                           dateo = ladate
                        else
                           call decodate( ncid,tim1, dateo )
                        endif
                        call decodate( ncid,tim2, ccctime )
                        call difdatr( ccctime,dateo, tdelta )
                        if (nint( dt ) <  1) then
                           npas = nint( tdelta ) ; deet = 3600
                           if (tim2-tim1 <= 1000000._8)         then
                               hold = tim2-tim1 
                               CALL Diag_CONVIP_plus
     +                         ( IP2,hold,kIND_HOURS,+2,NULS,.FALSE. )
                           else
                               IP2 = nint( tim2-tim1 )
                           end if
                           call            puthigh( IP2, 'IP2',ibuf )
                           svsm = 1 ; call puthigh( svsm,'IP3',ibuf )
                        else
                           loval = tim1
                           svsm = nint( tdelta / (dt/3600.0_8) )
                           npas = svsm-1 ; deet = nint( dt )
                           tdelta = dt / 3600.0_8 ; hival = tim2-tdelta
                           call incdatr( datei,dateo,tim1+tdelta )
                           dateo = datei
                        endif
                        call puthigh( dateo,'DATEO',ibuf )
                        call puthigh( npas, 'NPAS', ibuf )
                        call puthigh( deet, 'DEET', ibuf )
                        call puthigh( rkind,'RKIND',ibuf )
                        call puthir( hival, 'HIVAL',ibuf )
                        call puthir( loval, 'LOVAL',ibuf )
                        call putsamplz( svsm, ibuf )
                     else
                        call decodate( ncid,tim1,ccctime )
                        if (nint( dt ) >=  1) then
                           dateo = ladate
                           call difdatr( ccctime,dateo, tdelta )
                           npas = nint( tdelta / (dt/3600.0_8) )
                           deet = nint( dt )
                           call puthigh( dateo,'DATEO',ibuf )
                           call puthigh( npas, 'NPAS', ibuf )
                           call puthigh( deet, 'DEET', ibuf )
                        endif
                     endif
                  else
                     call decodate( ncid,tim1, ccctime )
                  endif
               else
                  ccctime=0
               endif

               ! Si _FillValue=NaN, ne pas arreter dans InfNaN
               holdl = noabort ; if (fill_cdf_nan) noabort = .true.

               if (timeid == -1) then
                  call get_vard(ncid,id,var(nn)%type,len,variable(1,nn))
               else
                  start(timeid)=itime
                  call get_varda(ncid,id,var(nn)%type,start,count,len,
     .                                                   variable(1,nn))
               endif

               noabort = holdl ! Restaurer noabort a sa valeur initiale

*              Traitement des donnees :

               do kk=nlev,1,-1                         ! Commence en haut
                  
                  indice=0
                  k=nlev-kk+1

                  fill_ccc_oui=.false.

                  if (transpose_xy) then
                     indice1 = dim1*dim2*(kk-1)
                     do i=1,dim1
                        do j=1,dim2
                           indice1 = indice1+1
                           indice2 = (j-1)*dim1+i
                           dval(indice2) = variable(indice1,nn)
                        end do
                     end do
                     indice1 = dim1*dim2*(kk-1)
                     variable(indice1+1:indice1+indice2,nn) =
     .               dval(1:indice2)
                  endif

                  call combline6( variable(1,nn),dval,indice,
     .                      kk, dim1,dim2,nlev, scale,offset,
     .                      var(nn)%mult,var(nn)%add,
     .                      fill_ccc_def,fill_cdf, invj,iii,
     .                      fill_all,fill_cdf_nan, opack3 )
                           
                  if (fill_all) then
                     if (rpn_info) then
                        call fill_high( ibuf,
     .                                  GRTYP,IG1,IG2,IG3,IG4,etiket )
                        call setlab(ibuf,type,ccctime,cccname,ilevel(k),
     .                              dim1+iii,dim2,nhem,opack1)
                        call prtlab2(' Passer (FILL-ALL) ...',
     .                                 ibuf )
                     endif
                     cycle
                  endif

                  if (fill_ccc_oui) fill_count = fill_count+1

                  if ( fill_ccc_oui
     .           .and. cdf2_mode == 'cdf2rpn'
     .           .and. opack1    >     -32  )   then
                     if (fill_message(id)) then
                        fill_message(id) = .false.
CCC                     write(6,6070) cccname,cloche
                     endif
CCC                  opack2 = -32
                     opack2 = opack1
                  else
                     opack2 = opack1
                  endif

                  if (opack2 <= -32) then
                     datyp=IEEE
                  elseif (opack2 == -16) then
                     datyp=TURBO16
                  else
                     datyp=STD
                  endif

                  call puthigh( datyp,'DATYP',ibuf )

*                 Ecriture dans le fichier de sortie :

C                 if (list(id)%name == 'ap' .or.
C    .                list(id)%name == 'b' )
C    .                write(6,'(3A4,A)') Dtyp,type,ibuf(1),GRTYP

                  call fill_high( ibuf, GRTYP,IG1,IG2,IG3,IG4,etiket )
                  call setlab(ibuf,type,ccctime,cccname,ilevel(k),
     .                       dim1+iii,dim2,nhem,opack2)

                  call fool_optimizer( Dtyp,type,ibuf(1) )

C                 if (list(id)%name == 'ap' .or.
C    .                list(id)%name == 'b' )
C    .                write(6,'(3A4,A)') dtyp,type,ibuf(1),GRTYP

                  call putfld2(funit,dval,ibuf,maxpk)
                  if (((cdf2_mode.eq.'cdf2rpn'  .and. RPN_INFO) .or.
     .                  cdf2_mode.eq.'cdf2ccc')                
     .           .and. itime.eq.dim(timedid)%len)
     .                                      call prtlab( ibuf )
               end do                                    !level
            end do                                      !time
         end if
      end do TRAITER_VARIABLES

      if (fill_ccc_def) then
        if (fill_count.ne.0) then     !SI on a assigne des valeurs de remplissage 
          write(6,1003) fill_ccc
        else
          write(6,1013)
        endif
      endif

*-----------------------------------------------------------------------
 1003 format(//
     .    '**VALEURS DE REMPLISSAGE DANS SORTIE CCC/RPN ASSIGNEES A : ',
     .      1pe12.5/)
 1013 format(//
     .   '**AUCUNE VALEUR DE REMPLISSAGE ASSIGNEE DANS SORTIE CCC/RPN'/)

 6000 format(//' VARIABLES :')
 6001 format(/" A l'appel, definir :   -cle_nhem avec valeur de 1",
     .        ' pour HEMIS NORD ou 2 pour HEMIS SUD car grille',
     .        ' polaire-stereographique et direction=cccma'/)
 6002 format(/' Inconsistance des dimensions'/
     .        ' ntime est ',I8,' et Dim#1*Dim#2*Dim#3 est ',I8,
     .        ' mais la taille totale est ',I8)
 6004 format(/' Ordre des dimensions non supporte'/
     .        ' DimXid, DimYid, DimZid = ',3I2,' plutot que 1 2 3'/)
 6012 format( ' dim#',A,' = ',I8/)
 6022 format( ' dim#',A,' est indefinie'/)
 6003 format( ' dim ',A,' introuvable'/)

 6020 format(/)
 6050 format(/a1,'*** La variable "_FillValue" n''existe pas dans le'
     .     / '    fichier netcdf: on remet fill_ccc_def a .false.')
 6060 format(/a1,'*** Le fichier netcdf contient des grilles lat/lon',
     .      ' a resolution uniforme')
 6061 format(/a1,'*** Le fichier netcdf contient des grilles qui ne',
     .      ' sont soit ni gaussiennes, soit ni a resolution uniforme')
 6062 format(/' Grilles regionale non-supportees: les latitudes sont ',
     .        ' inversees.'/' La variable logique INVJ est re-definie',
     .        ' a ',L1,A)
 6063 format(/' Pour que les latitudes soient croissantes,',
     .        ' invj doit etre re-defini a ',L1,A/)
 6064 format(/a1,'*** Le fichier netcdf contient des grilles qui',
     .      ' sont (de type Z) non tournee et etirees')
 6065 format(/a1,'*** Le fichier netcdf contient des grilles qui',
     .      ' sont (de type Y)')
 6070 format(/' Traitement des valeurs manquantes actif pour variable ',
     .        a,', npack redefini a 32 bits.',a1)
 6080 format(/' Variable sans reperes geographiques ',A/)

 6100 format(/' Parametres de grilles polaire-stereographiques...'/
     .        ' pi,pj,d60,dgrw = ',4e12.5/)
 6101 format(/' Parametres de grilles lat-lon regionales...'/
     .        ' olat,olon,dlat,dlon = ',4e12.5/)
*-----------------------------------------------------------------------
      end
      subroutine def_lon_lat( cid, ocoord,dcoord, nijla,
     .                        dcons,cincreas,glb,bgrd,gauss )

      implicit none

      include 'cdf2ccc.h'
      include 'infomem.h'
      include 'varmem.h'
      include 'dimmem.h'

*     Bernard Dugas  mai 2007
*
*     Determiner les valeurs de ocoord=min( coord ) et de
*     ocoord=delta( coord ) de la coordonnee 1D associee a
*     #cid et de certaines autres proprietees de cette
*     coordonnee telle que:

*      1) Est-elle a resolution constante      (dcons); 
*      2) Valeurs croissantes ou decroissantes (cincreas);
*      3) Couverture globale ou locale         (glb);
*      4) Grille de type A ou B                (bgrd);
*      5) Latitudes (donc cid=yid) gaussiennes (gauss).


*     Arguments

      integer, intent(in)  :: cid
      real,    intent(out) :: ocoord,dcoord
      logical, intent(out) :: dcons,cincreas,bgrd,gauss
      integer, intent(out) :: nijla,glb

*     Variables locales

      
      integer i,j,nlen,ilath
      real*8, dimension(:), allocatable :: sl,wl,cl,rad,wossl,cvaleur
      real*8  dcoordmax,dcoordmin, deltcoord, sum,sum2, hold,
     .        coordmax,coordmin, petit,rad2deg

*-----------------------------------------------------------------------

      rad2deg  = 90./asin(1d0)
      petit    = 0.0001 ! = 0.001 (autre valeur possible)

      nijla    =  2

      cincreas = .true.
      glb      =  -1
      bgrd     = .false.
      gauss    = .false.

      nlen     = dim(coord(cid)%dimid(1))%len

      allocate( cvaleur(nlen) )
      cvaleur(1:nlen) =  dcoordonne(1:nlen,cid)
      
      if (cid == xid) then
         ! longitudes monotones/positives
         do i=1,nlen
            if (cvaleur(i) < -0.000001)
     .          cvaleur(i) = cvaleur(i)+360.
         enddo
         do i=2,nlen
            if (cvaleur(i) < cvaleur(i-1))
     .          cvaleur(i) = cvaleur(i)+360.
         enddo
      else if (cid == yid .and. invj) then
         do i=1,nlen/2
            hold              = cvaleur(i)
            cvaleur(i)        = cvaleur(nlen-i+1)
            cvaleur(nlen-i+1) = hold
         enddo
      endif

      deltcoord = abs( cvaleur(2)-cvaleur(1) )
      dcoordmax = deltcoord ; dcoordmin = dcoordmax
      sum = deltcoord ; sum2 = sum*sum

      coordmax  = max( cvaleur(1),cvaleur(2) )
      coordmin  = min( cvaleur(1),cvaleur(2) )

      do i=2,nlen-1

         deltcoord = abs( cvaleur(i+1)-cvaleur(i) )

         sum = sum+deltcoord
         sum2 = sum2 + (deltcoord*deltcoord)

         if (deltcoord.lt.(dcoordmin-petit)) then
            nijla = 2
         else if (deltcoord.lt.(dcoordmin+petit)) then
            nijla = nijla+1
         endif

         dcoordmax = max( deltcoord,dcoordmax )
         dcoordmin = min( deltcoord,dcoordmin )

         coordmax  = max( cvaleur(i+1),coordmax )
         coordmin  = min( cvaleur(i+1),coordmin )

      enddo

      if (nlen > 1) then
         sum = sum/(nlen-1) ; sum2 = sum2/(nlen-1)
      endif

      dcoord = sum
      sum2 = sqrt( max( 0.0_8,sum2-(sum*sum) ) )

      if (dcoordmax-dcoordmin.lt.petit*(dcoordmax+dcoordmin)) then
         nijla = nlen
!!!      dcoord = 0.5*(dcoordmax+dcoordmin)
         dcons = .true.
      else
!!!      dcoord = dcoordmin
         dcons = .false.
      endif

      ocoord = coordmin

      if (cvaleur(1) .gt. cvaleur(2)) cincreas = .false.

      if (dcons .and. cid.eq.xid) then

         if (abs( coordmin ) .lt. petit) then

*           repetition du meridien de Greenwich (grille B) ?

            if (abs( coordmax - 360. ) .lt. petit) bgrd = .true.

*           verifier si la coordonnee va de 0 a 360

            if (bgrd .or. abs( coordmax + dcoord - 360. ) .lt. petit)
     .          glb = 0

         endif

      elseif (cid.eq.yid) then

         if (dcons) then

*     verifier que les extremites du domaine contiennent ou bien
*     1) les deux poles, ou bien 2) un des poles et une latitude
*     "proche" de l'equateur dans la meme hemisphere que ce pole.

*     Si l'une de ces conditions est verifiee, il s'agit peut-etre
*     d'une grille de type A ou B.

            if ((abs( coordmin          +90. ).lt.petit  .and.
     .           abs( coordmax          -90. ).lt.petit)        .or.
     .          (abs( coordmin-dcoord/2.+90. ).lt.petit  .and.
     .           abs( coordmax+dcoord/2.-90. ).lt.petit))       then
               glb = 0
            elseif
     .         ((abs( coordmin               ).lt.petit  .and.
     .           abs( coordmax          -90. ).lt.petit)        .or.
     .          (abs( coordmin-dcoord/2.     ).lt.petit  .and.
     .           abs( coordmax+dcoord/2.+90. ).lt.petit))       then
               glb = 1
            elseif
     .         ((abs( coordmin          +90. ).lt.petit  .and.
     .           abs( coordmax               ).lt.petit)        .or.
     .          (abs( coordmin-dcoord/2.+90. ).lt.petit  .and.
     .           abs( coordmax+dcoord/2.     ).lt.petit))       then
               glb = 2
            endif

            if ((glb.eq.0                                .and.
     .           abs( coordmin          +90. ).lt.petit  .and.
     .           abs( coordmax          -90. ).lt.petit)        .or.
     .          (glb.eq.1                                .and.
     .           abs( coordmin               ).lt.petit  .and.
     .           abs( coordmax          -90. ).lt.petit)        .or.
     .          (glb.eq.2                                .and.
     .           abs( coordmin          +90. ).lt.petit  .and.
     .           abs( coordmax               ).lt.petit))
     .
     .         bgrd = .true.

         else if (nlen.gt.1) then

            allocate( sl(nlen),wl(nlen),cl(nlen),rad(nlen),wossl(nlen) )

*           verifier pour des grilles gausiennes.

*           commencer par le cas d'une grille globale

            ilath = nlen/2

            call gaussg( ilath,sl,wl,cl,rad,wossl )
            call trigl2( ilath,sl,wl,cl,rad,wossl, 0 )

            rad = rad * rad2deg

            sum = 0

            if (cincreas) then
               do i=1,nlen
                  sum = sum + (cvaleur(i)-rad(i))**2
               enddo
            else
               do i=1,nlen
                  sum = sum + (cvaleur(nlen-i+1)-rad(i))**2
               enddo
            endif

            sum = sqrt( sum/nlen )

            if (sum.lt.petit) then
               glb = 0
               goto 100
            endif

*           grille hemisphere nord ?

            ilath = nlen

            call gaussg( ilath,sl,wl,cl,rad,wossl )
            call trigl2( ilath,sl,wl,cl,rad,wossl, 1 )

            rad = rad * rad2deg

            sum = 0

            if (cincreas) then
               do i=1,nlen
                  sum = sum + (cvaleur(i)-rad(i))**2
               enddo
            else
               do i=1,nlen
                  sum = sum + (cvaleur(nlen-i+1)-rad(i))**2
               enddo
            endif

            sum = sqrt( sum/nlen )

            if (sum.lt.petit) then
               glb = 1
               goto 100
            endif

*           grille hemisphere sud ?

            ilath = nlen

            call gaussg( ilath,sl,wl,cl,rad,wossl )
            call trigl2( ilath,sl,wl,cl,rad,wossl, 2 )

            rad = rad * rad2deg

            sum = 0

            if (cincreas) then
               do i=1,nlen
                  sum = sum + (cvaleur(i)-rad(i))**2
               enddo
            else
               do i=1,nlen
                  sum = sum + (cvaleur(nlen-i+1)-rad(i))**2
               enddo
            endif

            sum = sqrt( sum/nlen )

            if (sum.lt.petit) glb   = 2

 100        if (glb.ge. 0   ) gauss = .true.

            deallocate( sl,wl,cl,rad,wossl )

         endif

      endif

      deallocate( cvaleur )

      return

*-----------------------------------------------------------------------
      end
      subroutine fill_high( ibuf, GRTYP,IG1,IG2,IG3,IG4,etiket )

      implicit none

*     Bernard Dugas   mai 2007

      integer ibuf(*),IG1,IG2,IG3,IG4
      character*(*) GRTYP,etiket

      character label*12
*-----------------------------------------------------------------------

      label = etiket

      CALL puthic( 'NC'  ,'TYPVAR', IBUF )
      CALL puthic(  GRTYP ,'GRTYP', IBUF )

      CALL puthigh( IG1   ,'IG1'  , IBUF )
      CALL puthigh( IG2   ,'IG2'  , IBUF )
      CALL puthigh( IG3   ,'IG3'  , IBUF )
      CALL puthigh( IG4   ,'IG4'  , IBUF )

      CALL puthic(  label(1:04),'ETIK1', IBUF )
      CALL puthic(  label(5:08),'ETIK2', IBUF )
      CALL puthic(  label(9:12),'ETIK3', IBUF )

      return

*-----------------------------------------------------------------------
      end
      logical function checkres( initial,delta )

      implicit none

!     Bernard Dugas, jan 2012 - Checks for appropriate L-grid resolutions

      real     initial,delta

      checkres = .false.

      if (delta <= 0.0) return

      if ((nint(initial*100) *100  /= nint(initial*10000)  ) .or.
     .    (delta <  1.  .and.
     .     nint(delta  *1000)*1000 /= nint(delta  *1000000)) .or.
     .    (delta >= 1.  .and. delta <= 20. .and.
     .     nint(delta  *100) *100  /= nint(delta  *10000)  ) .or.
     .    (delta >  20. .and. delta <= 55. .and.
     .     nint(delta)       *1000 /= nint(delta  *1000)   )) return

      checkres = .true.

      return

!-----------------------------------------------------------------------
      end
      subroutine decode_single_level( ilevel,level_value )

      implicit none

!     Bernard Dugas, jul 2012 - Try to decode a single level attribut

      integer ilevel
      character*(*) level_value

      character fmt8*8,nombre*2
      character(len=128) clocal,debut,fin
      real      valeur
      integer   i,ib,vkind
!-----------------------------------------------------------------------

      clocal = level_value
      if (clocal == ' ') return

      ib    = index( clocal,' ')
      debut = clocal(1:ib-1)
      fin   = clocal(ib+1:128)

      write(nombre,'(i2.2)') min( max( ib-1, 10 ), 20 )
      fmt8='(E'//nombre//'.0)'
      read(debut,fmt8) valeur

      if (ib == len_trim( clocal )+1) then

        vkind=3 ! Type arbitraire
        call convpr( ilevel, valeur, vkind, +1 )

      else

         do i=1,128-ib
            if (fin(i:i) /= ' ') exit
         enddo
         fin = fin(i:128-ib)

         if (fin == 'm' ) vkind=0
         if (fin == 'sg') vkind=1
         if (fin == 'mb') vkind=2
         if (fin == 'ar') vkind=3
         if (fin == 'hy') vkind=5
         if (fin == 'th') vkind=6

         call convpr( ilevel, valeur, vkind, +1 )

      endif

      return

!-----------------------------------------------------------------------
      end
      logical function monotone( variable, ni,nj, index,upward )

      implicit none

      ! Bernard Dugas, apr 2016 - Check for monotonicity and its direction 

      integer :: index     ! Which indicy to check, first (1) or second (2) ?
      logical :: upward    ! If monotonic, is the variable increasing ?
      integer :: ni,nj
      real    :: variable( ni,nj )

      ! Local variables

      logical :: ok
      integer :: i,j
!-----------------------------------------------------------------------

      ok = .true.

      if (index == 1) then

         do j=1,nj
            do i=2,ni
               if (variable(i-1,j) >= variable(i,j)) then
                  ok = .false.
                  exit
               endif
            enddo
         enddo

      else if (index == 2) then

         do j=2,nj
            do i=1,ni
               if (variable(i,j-1) >= variable(i,j)) then
                  ok = .false.
                  exit
               endif
            enddo
         enddo

      else ! Unrecognized index value

         ok = .false. ;  upward = ok ; monotone = ok ; return

      endif

      if (ok) then  ! Monotonicaly increasing values

         upward = ok ;  monotone = ok

      else ! Otherwise, check for monotonicaly decreasing values

         upward = ok ; ok = .true.

         if (index == 1) then

            do j=1,nj
               do i=1,ni-1
                  if (variable(i,j) <= variable(i+1,j)) then
                     ok = .false.
                     exit
                  endif
               enddo
            enddo

         else

            do j=1,nj-1
               do i=1,ni
                  if (variable(i,j) <= variable(i,j+1)) then
                     ok = .false.
                     exit
                  endif
               enddo
            enddo

         endif

         monotone = ok
         
      endif

      return

!-----------------------------------------------------------------------
      end

