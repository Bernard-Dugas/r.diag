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
      subroutine get_coord2 (NCID)

      implicit none

      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'

      integer ncid

******
*
* AUTEUR Guy Bergeron       juin 2003
*
*     Identifie les coordonnees a partir de la du type derive list(id)
*     et defini le type derive coord(id) associe a chaque coordonnee
*     trouvee. Defini aussi les id (xid,yid,zid,tid) correspondant
*     a chaque coordonnee.
*
* REVISIONS
*
*  B.Dugas fevrier '16 :
*  - Reconnaitre long_name="model level number" pour la coordonnee levelist
*  B.Dugas mars '15 :
*  - Introduire la variante 'ht' de la coordonnee de hauteur en km.
*  B.Dugas aout '13 :
*  - Corriger le traitement des variables xcoord, ycoord et zcoord
*  - Corriger le traitement des "long_name" qui pourront etre associees
*    au type de coordonnee verticale "Arbitrary Levels"
*  B.Dugas aout '12 :
*  - Ajouter 'height' a la liste des coordonnees verticale reconnues
*  - Si hold_desc est non-nul et different de level_desc a la fin
*    de la routine, verifier que le nouveau level_desc a une valeur
*    acceptable/possible. Sinon, restaurer la valeur hold_desc
*  B.Dugas juillet '12 :
*  - Ajouter le support de 'Log Pressure Hybrid Levels' (VKIND=5002)
*  - ajouter 'sigma' dans la liste des coordonnees verticale reconnues
*  B.Dugas avril '12 :
*  - Si 'air_pressure' (sn) ou 'pressure' (ln), level_desc='Pressure Levels'
*  - Lorsque level_desc corresponds a une coordoonne qui suit le terrain,
*    coord(zid)%mult=1000. seulement si Is_CCC=T
*  B.Dugas decembre '11 :
*  - Ajouter 'hybrid height coordinate' (pour HadGEM2-ES)
*  B.Dugas mars '10 :
*  - Ajouter 'hybrid' comme descripteur possible de la coordonnee verticale
*  B.Dugas mai '09 :
*  - Coordonnees NetCDF "[xyzt]coord" specifiees en arguments ?
*  B.Dugas fev '09 :
*  - Coordonnees verticale en unite de Pa (et non hPa)
*  B.Dugas oct '08 :
*  - Toujours definir xdid,ydid,zdid == list(id)%dimid(1)
*  - Traiter les niveaux de type "Arbitrary Levels" (kind=3)
*  - Supporter la coordonnee verticale ayant pour nom "lev"
*  - Ajouter le support de coordonnees de niveau de sol et TOA
*  Bernard Dugas fevrier 2008 :
*  - S'assurer que les majuscules/minuscules dans le
*    nom des coordonnees ne soient plus significatives
*  - On cherche les attributs 'long_name' et 'standard_name'
*    des coordonnees afin de mieux les identifier
*  - La coordonnee verticale peut maintenant etre identifiee par
*    l'attribut 'nlevs' en plus des attributs 'z' et level'
*  - Effectuer une premiere definition de la variable 'level_desc'
*  - Ajouter le support des coordonnees verticales de hauteur
*    et de pression hybride
*  - Les attributs de coordonnees horizontales 'rlon' et 'rlat'
*    associes aux reperes geographiques tournes sont supportes
*  - S'assurer que les coordonnees ne sont definies qu'une seule
*    fois, eg. la presence de 'LON' et 'lon' dans un meme fichier
*    generera une erreur fatale
*  Anne Frigon aout 2006 : Ajoute test pour level_desc avec units "hPa" 
*                          car modifie units de "millibar" a hPa"
*                          pour level_desc='Pressure Levels' dans attribut_coord.f
*  Guy Bergeron Mai 2005 : Modification du test "IF" sur "lon" et "lat".
*
******

      logical  Is_CCC, ReConnu
      logical, save :: defx=.false.,defy=.false.,defz=.false.
      character(len=128) cfield,cfield2,hold_desc
      integer i,ii,nn,ln,sn,id,cid,nlen
      
      data ii,nn /max_attrs,max_attrs/   !empecher de defoncer la memoire
*-----------------------------------------------------------------------
      write(6,6000)

      Is_CCC = .false. ; if (ccc_pktyp(1:2) == 'PK') Is_CCC = .true.

      do id=1,nvars

*     Identifier les coordonnees :

      if (list(id)%ndim.eq.1.and.
     .                list(id)%name.eq.dim(list(id)%dimid(1))%name) then

         list(id)%var_ok=.false.

         call affecte_coord(ncoord,cid,list(id)%type,list(id)%ndim,
     .                                     list(id)%dimid,list(id)%name)
         coord(cid)%nattr=list(id)%nattr

         ii = -1 ; nn = -1 ; ln = -1 ; sn = -1

         call get_attribut(ncid,id,coord(cid)%nattr,coord(cid)%name)

         do i=1,coord(cid)%nattr
            if (attr(i)%name.eq.'axis'         ) ii=i
            if (attr(i)%name.eq.'units'        ) nn=i
            if (attr(i)%name.eq.'long_name'    ) ln=i
            if (attr(i)%name.eq.'standard_name') sn=i
         enddo

*        Identifier les id des coordonnees x,y,z,t (i.e. xid,yid,zid et tid) :

*                     =======================
         call up2low( coord(cid)%name, cfield )
*                     =======================

         if (cfield          == 'p'         .or. 
     .       cfield          == 'lev'       .or. 
     .       cfield          == 'plev'      .or. 
     .       cfield          == 'level'     .or.
     .       cfield          == 'ht'        .or. 
     .       cfield          == 'height'    .or. 
     .       cfield          == 'sigma'     .or. 
     .       cfield          == 'hybrid'    .or. 
     .       cfield          == 'levelist'  .or. 
     .       cfield          == 'nlevels'   .or. 
     .       coord(cid)%name ==  zcoord     .or. 
     .       attr(ii)%cvalue == 'z'         .or. 
     .                                      attr(ii)%cvalue.eq.'Z') then 
            zid=cid
            zdid=list(id)%dimid(1)

            if (defz) then
               write(6,6002) 'Z'
               call xit(' Get_Coord2 ',-2 )
            else
               defz = .true.
            endif

         else
     .   if (cfield         .eq.'x'         .or. 
     .       cfield         .eq.'x_2'       .or. 
     .       cfield         .eq.'lon'       .or. 
     .       cfield         .eq.'rlon'      .or. 
     .       cfield         .eq.'longitude' .or. 
     .       coord(cid)%name.eq. xcoord     .or. 
     .       attr(ii)%cvalue.eq.'x'         .or.                       
     .                                      attr(ii)%cvalue.eq.'X') then
            xid=cid
            xdid=list(id)%dimid(1)
            lon=trim(coord(cid)%name)

            if (defx) then
               write(6,6002) 'X'
               call xit(' Get_Coord2 ',-3 )
            else
               defx = .true.
            endif

         else
     .   if (cfield         .eq.'y'         .or. 
     .       cfield         .eq.'y_2'       .or. 
     .       cfield         .eq.'lat'       .or.
     .       cfield         .eq.'rlat'      .or. 
     .       cfield         .eq.'latitude'  .or. 
     .       coord(cid)%name.eq. ycoord     .or. 
     .       attr(ii)%cvalue.eq.'y'         .or.
     .                                      attr(ii)%cvalue.eq.'Y') then
            yid=cid
            ydid=list(id)%dimid(1)
            lat=trim(coord(cid)%name)

            if (defy) then
               write(6,6002) 'Y'
               call xit(' Get_Coord2 ',-4 )
            else
               defy = .true.
            endif

         else
     .   if (cfield         .eq.'time'      .or. 
     .       cfield         .eq.'t'         .or. 
     .       cfield         .eq. tcoord     .or. 
     .       attr(ii)%cvalue.eq.'t'         .or.
     .                                      attr(ii)%cvalue.eq.'T') then
            timedid=coord(cid)%dimid(1)
            tid=cid

            if (.not.no_time) then
               write(6,6002) 'Z'
               call xit(' Get_Coord2 ',-5 )
            endif

            no_time =.false.           

         endif

*NOTA: Definir "mult" et "add" de sorte que level=value*mult+add
*      soit egale a la sortie de "lvdcode". La consistance de l'encodage
*      et du decodage est assuree dans "def_level"
*
*      Nous ne connaisons pas encore level_desc. Nous pouvons
*      definir cette variable selon l'attribut 'standard_name'
*      ou 'long_name' de la coordonnee trouvee dans le fichier.
*      Si l'attribut non standart 'level_desc' est egalement
*      present, il sera relu plus tard par la routine
*      get_attribut.

         if(cid.eq.zid) then

*           Premiere definition de level_desc

            hold_desc = level_desc ; level_desc = ' ' ; ReConnu = .true.

            if (sn .ne. -1 ) then

*              Utilisons l'attribut 'standard_name'

               call clean_char( attr(sn)%cvalue,cfield,nlen )
               call up2low( cfield, cfield )

               if     (cfield == 'atmosphere_pressure_coordinate' .or.
     .                 cfield == 'air_pressure'                  ) then
                  level_desc = 'Pressure Levels'

               elseif (cfield == 'atmosphere_sigma_coordinate'   ) then
                  level_desc = 'Sigma Levels'

               elseif (cfield == 
     .              'atmosphere_hybrid_sigma_pressure_coordinate') then
                  level_desc = 'Hybrid Levels'

               elseif (cfield == 
     .          'atmosphere_hybrid_sigma_log_pressure_coordinate') then
                  level_desc = 'Log Pressure Hybrid Levels'

               elseif (cfield == 
     .              'atmosphere_hybrid_height_coordinate'        ) then
                  level_desc = 'Hybrid Height'

               elseif (cfield == 'height'                        ) then
                  level_desc = 'Height'

               elseif (cfield == 'gal-chen levels'               ) then
                  level_desc = 'Gal-Chen Levels'

               elseif (cfield == 'layers'           .or.
     .                 cfield == 'arbitrary_levels' .or.
     .                 cfield == 'top_of_atmosphere'.or.
     .                 cfield == 'soil_layers'                   ) then
                  level_desc = 'Arbitrary Levels'

               else
                  level_desc = cfield ; Reconnu = .false.

               endif

            endif

            if (ln /= -1 .and.
     .         (sn == -1 .or. (sn /= -1 .and. .not.ReConnu))) then

*              Utilisons l'attribut 'long_name' si la coordonnee
*              n'a pas deja ete reconnue dans le 'if' precedent

               call clean_char( attr(ln)%cvalue,cfield,nlen )
               call up2low( cfield, cfield )

               if     (cfield == 'level'           .or.
     .                 cfield == 'pressure'        .or.
     .                 cfield == 'pressure level') then
!!!  .                 cfield == 'pressure_level') then
                  level_desc = 'Pressure Levels'

               elseif (cfield == 'sigma'         ) then
                  level_desc = 'Sigma Levels'

               elseif (cfield == 'hybrid'        ) then
                  level_desc = 'Hybrid Levels'

               elseif (cfield == 'logp hybrid'    ) then
                  level_desc = 'Log Pressue Hybrid Levels'

               elseif (cfield == 'hybrid height coordinate') then
                  level_desc = 'Hybrid Height'

               elseif (cfield == 'height'           .or.
     .               cfield == 'ht_mean_of_range') then
                  level_desc = 'Height'

               elseif (cfield == 'gal-chen'      ) then
                  level_desc = 'Gal-Chen Levels'

               elseif (cfield == 'layers'             .or.
     .                 cfield == 'arbitrary levels'   .or.
     .                 cfield == 'top of atmosphere'  .or.
     .                 cfield == 'model level number' .or.
     .                 cfield == 'soil layers'   )    then
                  level_desc = 'Arbitrary Levels'

               else
                  level_desc = cfield ; ReConnu = .false.

               endif

            endif

*           Definissons coord(zid)%add et coord(zid)%mult
*           selon la valeur de l'attribut 'units'. On peut
*           aussi redefinir level_desc si cela n'a pas ete
*           precedemment fait avec une valeur reconnue

            coord(zid)%add =0.0               
            coord(zid)%mult=1.0

            if (nn > 0) then

               call clean_char( attr(nn)%cvalue,cfield,nlen )
               call up2low( cfield, cfield )

               if (cfield(1:5).eq."sigma") then
                  if (level_desc == ' ' .or. .not.ReConnu)
     .                level_desc = 'Sigma Levels'
                  if (Is_CCC) coord(zid)%mult=1000.0

               else if (cfield.eq."1") then
                  if (level_desc == ' ' .or. .not.ReConnu)
     .                level_desc = 'Hybrid Levels'
                  if (Is_CCC) coord(zid)%mult=1000.0

               else if (cfield.eq."hybrid_sigma_pressure") then
                  if (level_desc == ' ' .or. .not.ReConnu)
     .                level_desc = 'Hybrid Levels'
                  if (Is_CCC) coord(zid)%mult=1000.0
                  do i=1,coord(cid)%nattr
                     if (attr(i)%name.eq.'top_pressure')
     .                  Hyb_pt=attr(i)%rvalue(1)
                     if (attr(i)%name.eq.'reference_pressure')
     .                  Hyb_pref=attr(i)%rvalue(1)
                     if (attr(i)%name.eq.'exponent')
     .                  Hyb_r=attr(i)%rvalue(1)
                  enddo
                  if (cdf2_mode.eq.'cdf2rpn')
     .               call setpt( Hyb_pt,Hyb_pref,Hyb_r )

               else if (cfield.eq."hybrid_sigma_log_pressure") then
                  if (level_desc == ' ' .or. .not.ReConnu)
     .                level_desc = 'Log Pressure Hybrid Levels'
                  if (Is_CCC) coord(zid)%mult=1000.0
                  do i=1,coord(cid)%nattr
                     if (attr(i)%name.eq.'top_pressure')
     .                  Hyb_pt=attr(i)%rvalue(1)
                     if (attr(i)%name.eq.'reference_pressure')
     .                  Hyb_pref=attr(i)%rvalue(1)
                     if (attr(i)%name.eq.'exponent1')
     .                  Hyb_r=attr(i)%rvalue(1)
                     if (attr(i)%name.eq.'exponent2')
     .                  Hyb_r2=attr(i)%rvalue(1)
                  enddo
               else if (cfield     .eq."hpa"       .or.
     .                  cfield(1:4).eq."mbar"      .or.
     .                  cfield(1:8).eq."millibar") then
                  if (level_desc == ' ' .or. .not.ReConnu)
     .                level_desc = 'Pressure Levels'
                  coord(zid)%mult=1.0

               else if (cfield.eq."pa") then
                  if (level_desc == ' ' .or. .not.ReConnu)
     .                level_desc = 'Pressure Levels'
                  coord(zid)%mult=0.01

               else if (cfield ==  "m"   .or.
     .                  cfield == "km" ) then
                  if (level_desc == ' ' .or. .not.ReConnu)
     .                level_desc = 'Height'
                  coord(zid)%mult=1.0
                  if (cfield == "km" .or. 
     .               (Is_CCC .and. level_desc == 'Gal-Chen Levels'))
     .               coord(zid)%mult=1000.

               else if (cfield == "level"      .or.
     .                  cfield == "ordinal"    .or.
     .                  cfield == "arbitrary") then
                  if (level_desc == ' ' .or. .not.ReConnu)
     .                level_desc = 'Arbitrary Levels'
                  coord(zid)%mult=1.0

               elseif (level_desc /= 'Arbitrary Levels'
     .          .and.   hold_desc == ' ') then
                  write(6,6001)trim( attr(nn)%cvalue ),nlen,
     .            (ichar( attr(nn)%cvalue(i:i) ),i=1,nlen)
                  write(6,6101) "-lev 'Arbitrary Levels'"
                  call                            xit('get_coord2',-1)

               endif

            endif

            if (level_desc == ' ') then
               level_desc = hold_desc
            else if (level_desc /= hold_desc .and.
     .               hold_desc  /= ' '      ) then
               call up2low( level_desc,cfield )
               do i=1,nlvl
                  call up2low( possible%level(i),cfield2 )
                  if (cfield == cfield2) exit
               enddo
               if (i > nlvl) level_desc = hold_desc
            endif

         endif

      endif
      enddo
         
*-----------------------------------------------------------------------
 6000 format(/' COORDONNEES :')
 6001 format(/' GET_COORD2 : ',
     .       ' Unites de la coordonnee verticale non reconnue = ',a/
     .       ' Le codage des ',I3,' caracteres de cette valeur sont...'/
     .        15x,12(:,I4,:,I4,:,I4,:,I4,:,I4,:,I4,:,I4,:,I4,:,I4))
 6101 format(/" Il est possible d'eviter cette condition en utilisant"/
     .        ' "',a,'" sur la ligne de commande.'/)
 6002 format(/' GET_COORD2 : ',
     .        " Plus d'une coordonnee en ",a/
     .  ' S.V.P. Modifiez votre fichier source avec NCDUMP/NCGEN.'/)
 6003 format(/' level_desc                = ',A,1x,A/)
*-----------------------------------------------------------------------
      end
      block data get_coord2_data

      include 'cdf2ccc.h'
      include 'infomem.h'

      data  no_time /.true./

      end block data get_coord2_data
