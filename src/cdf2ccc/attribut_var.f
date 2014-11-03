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
      subroutine attribut_var(ID,nbr,MOT_CLE)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h' 
      include 'varmem.h'

      integer id,nbr
      character*(*) mot_cle

******
*
*AUTEUR Guy Bergeron   juin  2003
*
*     Definir les attributs associes a la variables var(id).name.
*     
*REVISIONS
*
*     Bernard Dugas juin 2013 :
*     - Definir les attributs "valid_range" et "_FillValue" lorsque
*       la variable fill_ccc_def est Vrai (donc lorsque le traitement
*       des valeurs manquantes a ete active). Ces attributs sont du
*       meme type que les variables auquelles ils sont associes.
*       Ceci implique en particulier que "_FillValue" pourrait
*       etre modifie si sa valeur est trop grande. 
*     Bernard Dugas aout 2012 :
*     - Tenir compte de la variable cell_method pour definir
*       le type de traitement temporel qu'on subit les donnees
*     Bernard Dugas juin et juillet2012 :
*     - Definir l'attribut 'actual_level' plutot que
*       'level_desc' lorsque le mode single_level est actif
*     - Ajouter le traitement des options single_level
*       et single_time via l'argument MOT_CLE
*     Bernard Dugas mai 2012 :
*     - Definir l'attribut cell_methods
*     Bernard Dugas fevrier 2009 :
*     - Modification a l'appel de affecte_attr
*     - nbits depends maintenant de var(id)%type
*     Bernard Dugas ete 2007 :
*     - Ajouter le support des grille a repaires geographiques tournes
*     - Remplacer l'attribut "grid_mapping_name" par "grid_mapping" et
*       remplacer l'attribut "coordinate" par "coordinates"
*         
******

      integer*1     i1dummy,i1range(max_len)
      integer*2     i2dummy,i2range(max_len)
      integer       idummy,irange(max_len)
      character*128 dummy
      real*4        rdum,srange(max_len)
      real*8        ddummy,drange(max_len)

      integer       nlen,nt,validi
      integer*2     valide(max_len)
      character*128 lename,string
      real*8        scale,offset,dim_range(max_len)
      integer       nbits, i,ii,iii, largest,addid,scalid

      real*8        bigfloat,bigint2,bigint1

      LOGICAL           rpn_info,rpn_debug,noabort
      COMMON  /ZZVERBO/ rpn_info
      COMMON  /ZZDEBUG/          rpn_debug
      COMMON  /ZZABORT/                    noabort

*-----------------------------------------------------------------------
      nbr=0 ; valide=0

      bigint1 = huge( i1dummy )
      bigint2 = huge( i2dummy )
      bigfloat = huge( rdum )

*1    Prendre les attributs du fichier attribut_netcdf.dat

      nt=len(var(id)%name)
      call def_name (var(id)%name,nt,'',lename,nlen)             ! Enlever les blancs
      lename=lename(1:nlen)

      if (mot_cle.ne.'projection') then

         call get_attr(lename,nbr,var(id)%mult,var(id)%add,var(id)%name)

         if (mot_cle.eq.'char') return
         
      else 

         if(var(id)%name.eq.'polar_stereographic') then

            call def_attr_polar_stereographic(nbr)

         elseif(var(id)%name.eq.'rotated_latitude_longitude'  .or.
     .          var(id)%name.eq.'rotated_pole'              ) then

            call def_attr_rotated_lat_lon(nbr)

         else
            write(6,6001) trim(var(id)%name)
            call                    xit('attribut',-1)
         endif
         return
      endif

      call conv_unit(range,var(id)%mult,var(id)%add,2)    ! convertir en unites netcdf

*2    (peut-etre) definir "valid_range" et "_fill_value" :

      if (.not.spec) then

         scale=1.0 ; offset=0.0 ; validi = 0

         if (var(id)%type == nf_byte .or. var(id)%type == nf_short) then

            if (var(id)%type == nf_byte)  nbits = 8
            if (var(id)%type == nf_short) nbits = 16

                                        ! directement inspire
            largest=ibits(-1,0,nbits)/2 ! de paccrin avec nbits=16
                                        ! (Laprise et Fortin 1984)
            offset=range(1)
            scale=(range(2)-range(1))/dble(largest)

            do i=1,nbr
               if(attr(i)%name.eq."unpacked_valid_range") then 
                  validi=i
                  call valide_range(i,valide,scale,offset)
                  call affecte_attr(nbr,nf_short,'valid_range',2,
     .                                      dummy,i1dummy,valide,
     .                                        idummy,rdum,ddummy)
               endif
            enddo 

         else if (var(id)%type == nf_float) then

            srange(1:2) = range(1:2)

            do i=1,nbr
               if(attr(i)%name.eq."unpacked_valid_range") then 
                  validi=i
                  call affecte_attr(nbr,nf_float,'valid_range',
     .                       attr(i)%len,dummy,i1dummy,i2dummy,
     .                                    idummy,srange,ddummy)
               endif
            enddo 

         else if (var(id)%type == nf_double) then

            do i=1,nbr
               if(attr(i)%name.eq."unpacked_valid_range") then 
                  validi=i
                  call affecte_attr(nbr,nf_double,'valid_range',
     .                        attr(i)%len,dummy,i1dummy,i2dummy,
     .                                        idummy,rdum,range)
               endif
            enddo 

         endif

         if (validi == 0 .and. fill_ccc_def) then

            drange(1:2) = (range(1:2) - offset)/scale

            if (var(id)%type == nf_byte) then
               i1range(1:2) = int(drange(1:2))
               call affecte_attr(nbr,nf_byte,'valid_range',2,
     .                                 dummy,i1range,i2dummy,
     .                                    idummy,rdum,ddummy)
            else
     .      if (var(id)%type == nf_short) then
               i2range(1:2) = int(drange(1:2))
               call affecte_attr(nbr,nf_short,'valid_range',2,
     .                                  dummy,i1dummy,i2range,
     .                                     idummy,rdum,ddummy)
            else
     .      if (var(id)%type == nf_float) then
               srange(1:2) = drange(1:2)
               call affecte_attr(nbr,nf_float,'valid_range',2,
     .                                  dummy,i1dummy,i2range,
     .                                   idummy,srange,ddummy)
            else
     .      if (var(id)%type == nf_double) then
               call affecte_attr(nbr,nf_double,'valid_range',2,
     .                                   dummy,i1dummy,i2dummy,
     .                                      idummy,rdum,drange)
            endif

         endif

         if (fill_ccc_def) then ! definir "_FillValue"

            drange(1) = (fill_ccc - offset)/scale ; validi = 0

            do i=1,nbr
               if(attr(i)%name.eq."missing_value") then
                  if(attr(i)%type.eq.nf_double) then 
                     attr(i)%dvalue(1) = drange(1)
                  else
     .            if(attr(i)%type.eq.nf_float) then
                     if (bigfloat < abs( drange(1) )) then
                        attr(i)%rvalue(1) = nf_fill_float
                        if (rpn_info)
     .                     write(6,6100) "missing_value",
     .                     trim( var(id)%name ),nf_fill_float
                     else
                        attr(i)%rvalue(1) = drange(1)
                     endif
                  else
     .            if(attr(i)%type.eq.nf_short) then
                     if (bigint2 < abs( drange(1) )) then
                        attr(i)%i2value(1) = nf_fill_int2
                        if (rpn_info)
     .                     write(6,6101) "missing_value",
     .                     trim( var(id)%name ),nf_fill_int2
                     else
                        attr(i)%i2value(1) = drange(1)
                     endif
                  else
     .            if(attr(i)%type.eq.nf_byte) then
                     if (bigint1 < abs( drange(1) )) then
                        attr(i)%i1value(1) = nf_fill_int1
                        if (rpn_info)
     .                     write(6,6101) "missing_value",
     .                     trim( var(id)%name ),nf_fill_int1
                     else
                        attr(i)%i1value(1) = drange(1)
                     endif
                  endif
               endif
            enddo 

            do i=1,nbr
               if(attr(i)%name.eq."_FillValue") then 
                  validi = i
                  attr(i)%type = var(id)%type
                  if(attr(i)%type.eq.nf_double) then
                     attr(i)%dvalue(1) = drange(1)
                  else
     .            if(attr(i)%type.eq.nf_float) then
                     if (bigfloat < abs( drange(1) )) then
                        attr(i)%rvalue(1) = nf_fill_float
                        if (rpn_info)
     .                     write(6,6100) "_FillValue",
     .                     trim( var(id)%name ),nf_fill_float
                     else
                        attr(i)%rvalue(1) =  drange(1)
                     endif
                  else
     .            if(attr(i)%type.eq.nf_short) then
                     if (bigint2 < abs( drange(1) )) then
                        attr(i)%i2value(1) = nf_fill_int2
                        if (rpn_info)
     .                     write(6,6101) "_FillValue",
     .                     trim( var(id)%name ),nf_fill_int2
                     else
                        attr(i)%i2value(1) = drange(1)
                     endif
                  else
     .            if(attr(i)%type.eq.nf_byte) then
                     if (bigint1 < abs( drange(1) )) then
                        attr(i)%i1value(1) = nf_fill_int1
                        if (rpn_info)
     .                     write(6,6101) "_FillValue",
     .                     trim( var(id)%name ),nf_fill_int1
                     else
                        attr(i)%i1value(1) = drange(1)
                     endif
                  endif
               endif
            enddo 

            if (validi == 0) then

               if (var(id)%type == nf_double) then
                  call affecte_attr(nbr,nf_double,'_FillValue',1,
     .                                     dummy,i1dummy,i2dummy,
     .                                        idummy,rdum,drange)
               else
     .         if (var(id)%type == nf_float) then
                  if (bigfloat < abs( drange(1) )) then
                     srange(1) = nf_fill_float
                     if (rpn_info)
     .                  write(6,6100) "_FillValue",
     .                  trim( var(id)%name ),nf_fill_float
                  else
                     srange(1) = drange(1)
                  endif
                  call affecte_attr(nbr,nf_float,'_FillValue',1,
     .                                    dummy,i1dummy,i2dummy,
     .                                     idummy,srange,ddummy)
               else
     .         if (var(id)%type.eq.nf_short) then
                  if (bigint2 < abs( drange(1) )) then
                     i2range(1) = nf_fill_int2
                     if (rpn_info)
     .                  write(6,6101) "_FillValue",
     .                  trim( var(id)%name ),nf_fill_int2
                  else
                     i2range(1) = drange(1)
                  endif
                  call affecte_attr(nbr,nf_short,'_FillValue',1,
     .                                    dummy,i1dummy,i2range,
     .                                       idummy,rdum,ddummy)
               else
     .         if (var(id)%type.eq.nf_byte) then
                  if (bigint1 < abs( drange(1) )) then
                     i1range(1) = nf_fill_int1
                     if (rpn_info)
     .                  write(6,6101) "_FillValue",
     .                  trim( var(id)%name ),nf_fill_int1
                  else
                     i1range(1) = drange(1)
                  endif
                  call affecte_attr(nbr,nf_byte,'_FillValue',1,
     .                                   dummy,i1range,i2dummy,
     .                                      idummy,rdum,ddummy)

               endif
            endif

         endif

      endif


*3    Definir actual_range :

      if (mot_cle.eq.'auxiliary') then

         dim_range(1)=vlarge ; dim_range(2)=-vlarge

         nlen=1
         do i=1,var(id)%ndim
            if(var(id)%dimid(i).ne.timedid)
     .                               nlen=nlen*dim(var(id)%dimid(i))%len
         enddo

         call minmax(variable(1,id),dim_range(1),dim_range(2),nlen) !actual range

         call affecte_attr(nbr,nf_double,'actual_range',2,dummy,i1dummy,
     .                                    i2dummy,idummy,rdum,dim_range)
        return

      else if (validi > 0 .or. fill_ccc_def) then

         call affecte_attr(nbr,nf_double,'actual_range',2,dummy,i1dummy,
     .                                       i2dummy,idummy,rdum, range)

      end if

*4    Definir add_offset et scale_factor :
*
      if (.not.spec .and.
     .   (var(id)%type == nf_byte     .or.
     .    var(id)%type == nf_short) ) then

         call affecte_attr(nbr,nf_double,'add_offset',1,dummy,i1dummy,
     .                                     i2dummy,idummy,rdum,offset)

         call affecte_attr(nbr,nf_double,'scale_factor',1,dummy,i1dummy,
     .                                        i2dummy,idummy,rdum,scale)
      endif


*5    Definir coordinate :
      
      if (mot_cle.ne.'auxiliary') then

         if(project%name.eq.'rotated_latitude_longitude' .or.
     .      project%name.eq.'rotated_pole'               .or.
     .      project%name.eq.'polar_stereographic'      ) then

            string='lon lat'
            nlen=len_trim(string)
            call affecte_attr(nbr,nf_char,'coordinates',nlen,string,
     .                            i1dummy,i2dummy,idummy,rdum,ddummy)

            string=project%name
            nlen=len_trim(string)
            call affecte_attr(nbr,nf_char,'grid_mapping',nlen,
     .                        string,i1dummy,i2dummy,idummy,rdum,ddummy)

         else if(project%name(1:7).ne.'lon/lat'   .and. 
     .           project%name     .ne.'gaussian') then
            write(6,6002) trim(var(id)%name)
            call                                   xit('attribut',-2)
         endif
      endif


*6    Definir level_desc :

      if (index( mot_cle,'single_level' ) > 0) then
         ii=index( mot_cle,'single_level' )
         iii=index( mot_cle,';' )
         if(iii == 0) iii=len_trim( mot_cle )+1
         string=mot_cle(ii+13:iii-1)
         nlen=len_trim( string )
         call affecte_attr(nbr,nf_char,'actual_level',nlen,string,
     .                          i1dummy,i2dummy,idummy,rdum,ddummy)
      else
         nlen=len_trim(level_desc)
         call affecte_attr(nbr,nf_char,'level_desc',nlen,level_desc,
     .                          i1dummy,i2dummy,idummy,rdum,ddummy)
      endif

*7    Definir time_desc :
                   
      if (index( mot_cle,'single_time' ) > 0) then
         ii=index( mot_cle,'single_time' )
         string=mot_cle(ii+12:len_trim( mot_cle ))
         nlen=len_trim( string )
         call affecte_attr(nbr,nf_char,'time_desc',nlen,string,
     .                          i1dummy,i2dummy,idummy,rdum,ddummy)
      endif

*8    Definir grid_desc :

      string=project%name
      if (string(1:14).eq.'lon/lat global') string='lon/lat global'

      nlen=len_trim( string )
      call affecte_attr(nbr,nf_char,'grid_desc',nlen,string,
     .                             i1dummy,i2dummy,idummy,rdum,ddummy)

*9    Definir cell_methods:

      if (time_bnds_L) then
         if (cell_method == '?') then
            string="time: mean"
         else
            string=cell_method
         endif
      else
         string="time: point"
      endif

      nlen=len_trim( string ) 
      call affecte_attr(nbr,nf_char,'cell_methods',nlen,string,
     .                      i1dummy,i2dummy,idummy,rdum,ddummy)   

*-----------------------------------------------------------------------
 6001 format(' ATTRIBUT_VAR : Les attributs ne sont pas definis pour ',
     .                                                                a)
 6002 format(' ATTRIBUT_VAR : Il faut definir les attributs ',
     .                      '"grid_mapping" et "coordinates"',
     .                      ' pour la variable ', A)
 6100 format(" ATTRIBUT_VAR INFO : La valeur de l'attribut ",A,
     .       ' de la variable ',A,' est redefini a ',E20.10)
 6101 format(" ATTRIBUT_VAR INFO : La valeur de l'attribut ",A,
     .       ' de la variable ',A,' est redefini a ',I8)
*----------------------------------------------------------------------
      end

