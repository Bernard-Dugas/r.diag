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
      subroutine attribut_coord (ID,nbr)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'
      include 'varmem.h'
      include 'workmem.h'

      integer id,nbr
      
      
******
*
*AUTEUR Guy Bergeron   juin 2003
*
*     Definissons les attributs relatif a la variables de coordonnnee
*     coord(id).
*     
*REVISIONS
*
*  Bernard Dugas, septembre 2014 :
*  - Remplacer le '0.0' a la fin des unites temporelles par un '00'
*  Bernard Dugas, aout 2012 :
*  - Ajouter 'height' a la liste des coordonnees verticale reconnues
*  - Le nom par defaut de la coordonnee Z passe de 'level' a 'lev'
*  Bernard Dugas, juillet 2012 :
*  - Re-Definir dim_range pour la variable temporelle
*  - Ajouter le support de 'Log Pressure Hybrid Levels' (VKIND=5002)
*    lorsque les time_bnds sont actifs.
*  Bernard Dugas, juin 2012 :
*  - Ajouter la coordonnee verticale THETA
*  - Ajouter les formules (AP,B) pour VKIND=1,5,5001 et 5002
*  - Faire la distinction entre les calendriers 'gregorian'/
*    'standard' et 'proleptic_gregorian': On ne peut utiliser
*    les fonctions basees sur UDUNITS avec 'proleptic_gregorian'
*  Bernard Dugas, mai 2012 :
*  - Introduit le support des 'time_bnds'
*  - Commenter l'appel a get_att pour coord(id)%name = 'level'
*  - Ajouter le support de 'plev' et 'height' pour les noms
*    de la coordonnee verticale
*  Bernard Dugas, decembre 2011 :
*  - Toujour definir l'attribut calendar. Les valeurs
*    connues sont 'gregorian', '365_day' et '360_day'.
*  - Ajouter le support de la coordonnee verticale 'Hybrid Height'.
*  Bernard Dugas, juillet 2011 :
*  - Si leap=.F., definir l'attribut calendar='365_day'
*  Bernard Dugas, fevrier 2009 :
*  - Modification a l'appel de affecte_attr
*  Bernard Dugas, octobre 2008:
*  - Ajouter le support de coordonnees arbitraire, niveau de sol et TOA
*  Bernard Dugas mai,juillet 2008 :
*  - Supporter delta%type={s,m,h,d,M et y}
*  - Valeurs de l'attribut "axis" en majuscules (X,Y,Z et T)
*  Bernard Dugas hiver 2007 :
*  - Ajouter le support des coordonnes de pression hybride et de
*    hauteur. Cette derniere est distincte des options '10 m'/'2 m'
*  - Faire un "call xit" en cas d'erreur
*  Anne Frigon aout 2006 : Modifie units de "millibar" a hPa" 
*                          pour level_desc='Pressure Levels'
*  Guy Bergeron avril 2004 : Declaration de coordonne en REAL*8
*
******

      integer i,ni,clen
      integer year,month,day,hour,minute,sec
      real    xjour,xheure,xminute,xsec
      real*8  dim_range(2)

      integer*1    i1dummy
      integer*2    i2dummy
      integer      idummy
      real         rdummy
      real*8       ddummy

      character*128 dummy,string
      integer nlen
*-----------------------------------------------------------------------
      nbr=0

      clen=dim(coord(id)%dimid(1))%len

      do i=1,clen
         dval(i)=dcoordonne(i,id) ! dval est utilisee pour definir dim_range
      enddo
      
      if(coord(id)%name /= 'lev'   .and.
     .   coord(id)%name /= 'plev'  .and.
     .   coord(id)%name /= 'height')
     .        call get_attr(coord(id)%name,nbr,coord(id)%mult,
     .                                              coord(id)%add,dummy)

      if(coord(id)%name.eq.'time') then    ! Attributs complementaire de time

         nbr=nbr+1
         attr(nbr)%type=nf_char       ! character

         attr(nbr)%name='axis'
         attr(nbr)%cvalue='T'
         ni=len(attr(nbr)%cvalue)
         call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len) !longueur de cvalue

*     definir l'attribut 'calendar'

         nbr=nbr+1
         attr(nbr)%type=nf_char ! character
         attr(nbr)%name='calendar'

         if (cccvx) then ! Calendrier 360 jours
            attr(nbr)%cvalue='360_day'
         else if (.not.leap) then ! Calendrier non-gregorien
            attr(nbr)%cvalue='365_day'
         else if (noUD) then
            attr(nbr)%cvalue='proleptic_gregorian'
         else
            attr(nbr)%cvalue='gregorian'
         end if

         ni=len(attr(nbr)%cvalue)
         call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len) !longueur de cvalue

*     syntaxe de l'attribut "units" :  "????? since 1800-1-1 00:00:0.0"    


         call def_date2( ladate, year,month,day,hour,minute,sec,
     .                                                         'decode')

         nbr=nbr+1
         attr(nbr)%type=nf_char       ! character

         attr(nbr)%name='units'

         if      (delta%type.eq.'s') then
                  string='seconds since '              

         else if (delta%type.eq.'m') then
                  string='minutes since '              

         else if (delta%type.eq.'h') then
                  string='hours since '              

         else if (delta%type.eq.'d') then
                  string='days since '              
            
         else if (delta%type.eq.'M') then
                  string='months since '              

         else if (delta%type.eq.'y') then
                  string='years since '              

         endif
         
         if (sec == 0) then
            write(attr(nbr)%cvalue,7777) trim( string ),
     .                             year,month,day,hour,minute,sec
CCC                                year,month,day,hour,minute,0.0
         else
            write(attr(nbr)%cvalue,7778) trim( string ),
     .                             year,month,day,hour,minute,sec
         endif

         ni=len(attr(nbr)%cvalue)
         call get_string (attr(nbr)%cvalue,ni,'}',dummy,attr(nbr)%len)

         attr(nbr)%cvalue=attr(nbr)%cvalue(1:attr(nbr)%len)        ! enlever le "}"


*     syntaxe de l'attribut "delta_t" : "0000-00-00 06:00:00"


         year = 0 ; month  = 0 ; day = 0
         hour = 0 ; minute = 0 ; sec = 0


         if      (delta%type.eq.'s') then

            xminute= delta%dval          /60.
            xheure = xminute             /60.
            xjour  = xheure              /24.

            xheure = xheure -nint(xjour) *24.
            xminute= xminute-nint(xheure)*60.

            sec    = mod(nint(delta%dval),60)
            minute = nint( xminute )
            hour   = nint( xheure )
            day    = nint( xjour )

         else if (delta%type.eq.'m') then

            xheure = delta%dval          /60.
            xjour  = xheure              /24.

            xheure = xheure-nint(xjour)  *24.

            minute = mod(nint(delta%dval),60 )
            hour   = nint( xheure )
            day    = nint( xjour )

         else if (delta%type.eq.'h') then

            xjour  = delta%dval    /24.

            hour   = mod( nint( delta%dval ),24 )
            day    = nint( xjour )

         else if (delta%type.eq.'d') then

            day    = nint( delta%dval )

         else if (delta%type.eq.'M') then

            year   = nint( delta%dval     /12 )
            month  = nint( delta%dval-year*12 )

         else if (delta%type.eq.'y') then

            year   = nint( delta%dval ) 

         endif

         nbr=nbr+1
         attr(nbr)%type=nf_char        ! character

         attr(nbr)%name='delta_t'
         write(attr(nbr)%cvalue,8888)'',year,month,day,hour,minute,sec  

         i=index(attr(nbr)%cvalue,'{')
         ni=len(attr(nbr)%cvalue)
         call get_string (attr(nbr)%cvalue,ni,'}',dummy,
     .                                              attr(nbr)%len)
         attr(nbr)%cvalue=attr(nbr)%cvalue(1:attr(nbr)%len)

*     avg_period :

         if (time_bnds_L) then

            nbr = nbr+1
            attr(nbr)%type=nf_char       ! character
            attr(nbr)%name='bounds'
            attr(nbr)%cvalue='time_bnds'

            ni=len(attr(nbr)%cvalue)
            call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len) !longueur de cvalue

         endif
      
*     actual_range :

         if (time_bnds_L) then
            dim_range(1) = 0.5*(time_bnds(1,1)+time_bnds(2,1))
            dim_range(2) = 0.5*(time_bnds(1,clen)+time_bnds(2,clen))
         else
            dim_range(1)=dval(1)
            dim_range(2)=dval(clen)
         endif

      else if(coord(id)%name == 'lev'   .or.
     .        coord(id)%name == 'plev'  .or.
     .        coord(id)%name == 'height') then    ! Attributs complementaires de level

         nbr=nbr+1
         attr(nbr)%type=nf_char       ! character
         attr(nbr)%name='units'

         if (level_desc.eq.'Pressure Levels') then
            attr(nbr)%cvalue='Pa'

         else if (level_desc.eq.'Sigma Levels') then
            attr(nbr)%cvalue='sigma'

         else if (level_desc.eq.'Theta Levels') then
            attr(nbr)%cvalue='K'

         else if (level_desc.eq.'Hybrid Levels') then
            if (vkind /= 5002) then
               attr(nbr)%cvalue='hybrid_sigma_pressure'
            endif

         else if (level_desc.eq.'Log Pressure Hybrid Levels') then
            if (vkind == 5002) then
               attr(nbr)%cvalue='hybrid_sigma_log_pressure'
            endif

         else if (level_desc.eq.'Arbitrary Levels') then
            attr(nbr)%cvalue='arbitrary'

         else if (level_desc.eq.'Soil Layers') then
            attr(nbr)%cvalue='layers'

         else if (level_desc.eq.'Top of Atmosphere') then
            attr(nbr)%cvalue='TOA'

         else if (level_desc.eq.'Gal-Chen Levels' .or.
     .            level_desc.eq.'Hybrid Height'   .or.
     .            level_desc.eq.'Height'          .or.
     .            level_desc.eq.'10 m'            .or.
     .            level_desc.eq.'2 m'  ) then

            attr(nbr)%cvalue='m' 
         else
            print*,level_desc,'attribut_coord : Il faut definir units'
            call xit('attribut_coord',-1)
         endif

         attr(nbr)%len=len_trim(attr(nbr)%cvalue)

         nbr=nbr+1
         attr(nbr)%type=nf_char       ! character
         attr(nbr)%name='long_name'

         if      (level_desc.eq.'Pressure Levels') then
            attr(nbr)%cvalue='pressure level'

         else if (level_desc.eq.'Sigma Levels'   ) then
            attr(nbr)%cvalue='sigma'

         else if (level_desc.eq.'Theta Levels'   ) then
            attr(nbr)%cvalue='theta'

         else if (level_desc.eq.'Hybrid Levels'  ) then
            attr(nbr)%cvalue='hybrid'

         else if (level_desc == 'Log Pressure Hybrid Levels'  ) then
            attr(nbr)%cvalue='logp hybrid'

         else if (level_desc.eq.'Arbitrary Levels') then
            attr(nbr)%cvalue='arbitrary levels'

         else if (level_desc.eq.'Soil Layers') then
            attr(nbr)%cvalue='soil layers'

         else if (level_desc.eq.'Top of Atmosphere') then
            attr(nbr)%cvalue='top of atmosphere'

         else if (level_desc.eq.'Gal-Chen Levels') then
            attr(nbr)%cvalue='Gal-Chen'                

         else if (level_desc.eq.'Hybrid Height') then
            attr(nbr)%cvalue='hybrid height coordinate'

         else if (level_desc.eq.'Height'         ) then
            attr(nbr)%cvalue='height'                

         else
            attr(nbr)%cvalue=level_desc
         endif

         attr(nbr)%len=len_trim(attr(nbr)%cvalue)

         nbr=nbr+1
         attr(nbr)%type=nf_char       ! character
         attr(nbr)%name='standard_name'

         if      (level_desc.eq.'Pressure Levels') then
            attr(nbr)%cvalue='air_pressure'

         else if (level_desc.eq.'Sigma Levels'   ) then
            attr(nbr)%cvalue='atmosphere_sigma_coordinate'

         else if (level_desc.eq.'Theta Levels'   ) then
            attr(nbr)%cvalue='atmosphere_theta_coordinate'

         else if (level_desc.eq.'Hybrid Levels'  ) then
            if (vkind /= 5002)
     .        attr(nbr)%cvalue=
     .       'atmosphere_hybrid_sigma_pressure_coordinate'
               
         else if (level_desc == 'Log Pressure Hybrid Levels') then
            if (vkind == 5002) 
     .        attr(nbr)%cvalue=
     .       'atmosphere_hybrid_sigma_log_pressure_coordinate'

         else if (level_desc.eq.'Hybrid Height') then
            attr(nbr)%cvalue=
     .     'atmosphere_hybrid_height_coordinate'

         else if (level_desc.eq.'Height'         ) then
            attr(nbr)%cvalue='height'

         else if (level_desc.eq.'Arbitrary Levels') then
            attr(nbr)%cvalue='arbitrary_levels'

         else if (level_desc.eq.'Soil Layers') then
            attr(nbr)%cvalue='soil_layers'

         else if (level_desc.eq.'Top of Atmosphere') then
            attr(nbr)%cvalue='top_of_atmosphere'

         else
            attr(nbr)%cvalue=level_desc
         endif

         attr(nbr)%len=len_trim(attr(nbr)%cvalue)

         if (level_desc == 'Hybrid Levels') then

            if (vkind ==    1  .or.
     .          vkind ==    5  .or. 
     .          vkind == 5001) then

               nbr=nbr+1
               attr(nbr)%type=nf_char ! character

               attr(nbr)%name='formula'
               attr(nbr)%cvalue='p(n,k,j,i)=ap(k)+b(k)*ps(n,j,i)'
               ni=len_trim(attr(nbr)%cvalue)
               call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len)

               nbr=nbr+1
               attr(nbr)%type=nf_char ! character

               attr(nbr)%name='formula_terms'
               attr(nbr)%cvalue='a: ap b: b ps: ps'
               attr(nbr)%len=len_trim(attr(nbr)%cvalue)

CCC            ni=len_trim(attr(nbr)%cvalue)
CCC            call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len)

               call affecte_attr(nbr,nf_float,'top_pressure',
     .                  1,dummy,i1dummy,i2dummy,idummy,Hyb_pt,ddummy)
               call affecte_attr(nbr,nf_float,'reference_pressure',
     .                1,dummy,i1dummy,i2dummy,idummy,Hyb_pref,ddummy)
               call affecte_attr(nbr,nf_float,'exponent',
     .                   1,dummy,i1dummy,i2dummy,idummy,Hyb_r,ddummy)

            endif

         else if (level_desc == 'Log Pressure Hybrid Levels') then

            if (vkind == 5002) then

               nbr=nbr+1
               attr(nbr)%type=nf_char ! character

               attr(nbr)%name='formula'
               attr(nbr)%cvalue=
     .                'ln(p(n,k,j,i))=ln(ap(k))+b(k)*ln(ps(n,j,i)/pref)'
               ni=len_trim(attr(nbr)%cvalue)
               call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len)

               nbr=nbr+1
               attr(nbr)%type=nf_char ! character

               attr(nbr)%name='formula_terms'
               attr(nbr)%cvalue=
     .                     'ap: ap b: b ps: ps pref: reference_pressure'
               attr(nbr)%len=len_trim(attr(nbr)%cvalue)

CCC            ni=len_trim(attr(nbr)%cvalue)
CCC            call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len)

               call affecte_attr(nbr,nf_float,'top_pressure',
     .                  1,dummy,i1dummy,i2dummy,idummy,Hyb_pt,ddummy)
               call affecte_attr(nbr,nf_float,'reference_pressure',
     .                1,dummy,i1dummy,i2dummy,idummy,Hyb_pref,ddummy)
               call affecte_attr(nbr,nf_float,'exponent1',
     .                   1,dummy,i1dummy,i2dummy,idummy,Hyb_r,ddummy)
               call affecte_attr(nbr,nf_float,'exponent2',
     .                  1,dummy,i1dummy,i2dummy,idummy,Hyb_r2,ddummy)

            endif

         endif

         nbr=nbr+1
         attr(nbr)%type=nf_char       ! character

         attr(nbr)%name='axis'
         attr(nbr)%cvalue='Z'
         ni=len(attr(nbr)%cvalue)
         call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len) !longueur de cvalue

*
         if (level_desc.eq.'Sigma Levels'               .or.
     .       level_desc.eq.'Hybrid Levels'              .or.
     .       level_desc.eq.'Log Pressure Hybrid Levels' .or.
     .       level_desc.eq.'Pressure Levels'            .or.
     .       level_desc.eq.'Arbitrary Levels'           .or.
     .       level_desc.eq.'Soil Layers'              ) then

            nbr=nbr+1
            attr(nbr)%type=nf_char    ! character
            attr(nbr)%name='positive'
            attr(nbr)%cvalue='down'

         else

            nbr=nbr+1
            attr(nbr)%type=nf_char    ! character
            attr(nbr)%name='positive'
            attr(nbr)%cvalue='up'

         endif

         ni=len(attr(nbr)%cvalue)
         call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len) !longueur de value
*
         dim_range(1)=vlarge ; dim_range(2)=-vlarge
         call minmax(dval,dim_range(1),dim_range(2),
     .                                      dim(coord(id)%dimid(1))%len) !actual range

CCC      call get_attr(coord(id)%name,nbr,ddummy,ddummy,dummy) !bypass de mult et add?


      else if(coord(id)%name.eq.trim(lon)) then

*     Longitude (lon):

         nbr=nbr+1
         attr(nbr)%type=nf_char       ! character
         attr(nbr)%name='axis'
         attr(nbr)%cvalue='X'
         ni=len(attr(nbr)%cvalue)
         call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len) !longueur de cvalue

         dim_range(1)=vlarge ; dim_range(2)=-vlarge
         call minmax(dval,dim_range(1),dim_range(2),
     .                                      dim(coord(id)%dimid(1))%len) !actual range

      else if(coord(id)%name.eq.trim(lat)) then

*     Latitude (lat):

         nbr=nbr+1
         attr(nbr)%type=nf_char       ! character
         attr(nbr)%name='axis'
         attr(nbr)%cvalue='Y'
         ni=len(attr(nbr)%cvalue)
         call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len) !longueur de cvalue

         dim_range(1)=vlarge ; dim_range(2)=-vlarge
         call minmax(dval,dim_range(1),dim_range(2),
     .                                      dim(coord(id)%dimid(1))%len) !actual range
         
      else if(coord(id)%name.eq.'xc'    .or.
     .        coord(id)%name.eq.'rlon') then

*     Coordonne X pour grille polaire-stereograpique
*                   ou bien
*     Coordonne X pour grille avec rotation :

         nbr=nbr+1
         attr(nbr)%type=nf_char       ! character
         attr(nbr)%name='axis'
         attr(nbr)%cvalue='X'
         ni=len(attr(nbr)%cvalue)
         call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len) !longueur de cvalue

         dim_range(1)=vlarge ; dim_range(2)=-vlarge
         call minmax(dval,dim_range(1),dim_range(2),
     .                                      dim(coord(id)%dimid(1))%len) !actual range
         
      else if(coord(id)%name.eq.'yc'    .or.
     .        coord(id)%name.eq.'rlat') then

*     Coordonne Y pour grille polaire-stereograpique
*                   ou bien
*     Coordonne Y pour grille avec rotation :

         nbr=nbr+1
         attr(nbr)%type=nf_char       ! character
         attr(nbr)%name='axis'
         attr(nbr)%cvalue='Y'
         ni=len(attr(nbr)%cvalue)
         call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len) !longueur de cvalue

         dim_range(1)=vlarge ; dim_range(2)=-vlarge
         call minmax(dval,dim_range(1),dim_range(2),
     .                                      dim(coord(id)%dimid(1))%len) !actual range
         
      end if
      
*     Coefficients spectraux (a faire) :
      
      if(spec) then
         print *,'attribut_coord: ne fonctionne pas pour cas spectral'
         call xit('attribut_coord',-2)
      endif

*****     
            nbr=nbr+1
            attr(nbr)%type=nf_char ! float
            attr(nbr)%name='coordinate_defines'
            attr(nbr)%cvalue='point'
            ni=len(attr(nbr)%cvalue)
            call def_name(attr(nbr)%cvalue,ni,'',dummy,attr(nbr)%len) !longueur de cvalue

      nbr=nbr+1
      attr(nbr)%type=nf_double         ! double
      attr(nbr)%name='actual_range'
      attr(nbr)%len =2
      do i=1,2
         attr(nbr)%dvalue(i)=dim_range(i)
      end do
*-----------------------------------------------------------------------
      include 'format.h'
 9997 format('Hybrid model level Pt,P0=',f10.6,1x,f5.0,
     .      ' (Pa) and r=',f4.2)
 9998 format(a15,2(i3),a)           !attrc
 9999 format(a15,2(i3),2(f12.3))    !attrf

*-----------------------------------------------------------------------
      end

