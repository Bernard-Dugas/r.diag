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
      subroutine define_netcdf2( ncid,FUNIT,PHIS_UNIT,IBUF )

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'

      integer funit,phis_unit,ibuf(1)
      integer ncid

******
*
*AUTEUR Guy Bergeron         juin 2003
*
*     Traduction : CCCma -> netCDF 
*
*     A partir du fichier CCCma identifier les dimensions (i.e. dim(i)) et les 
*     valeurs des coordonnees (i.e. coord(i)). Creer le fichier netCDF. En mode
*     "define" dans le fichier netCDF, definir les dimensions, les coordonnees 
*     et la variables en plus d'associer les attributs aux coordonnees et aux 
*     variables. 
*
*REVISIONS
*
*
*  Bernard Dugas juin/juillet 2013 :
*  - Tenir compte de la variable meta_title -> attribut global 'title'
*  - Seulement enlever le mode "fill mode" lorsque la variable
*    fill_ccc_def est Faux (pas de traitement des valeurs manquantes)
*  Bernard Dugas aout 2012 :
*  - Meme si nlev=1, ne pas definir single_level
*    lorsque infvar(:)%unique_L est vrai
*  Bernard Dugas juin/juillet 2012 :
*  - La variable string passe a 512 caracteres et on corrige
*    le numero de l'attribut 'title' de la variable globale
*  - Tenir compte explicitement de 'Log Pressure Hybrid Levels'
*  - Ne pas activer le mode 'single_level' pour les fichiers CCC
*  - Allouer asser d'espace dans mot_cle pour tout definir (BuxFix)
*  - Transmettre le niveau lui-meme (en format caractere decode)
*    a attribut_var lorsque infvar(:)%len(zid)=1
*  - Transmettre la date elle-meme (en format caractere)
*    a attribut_var lorsque infvar(:)%len(tid)=1
*  - Definir les attributs associes aux vecteurs AP
*    et B pour les coordonnees de pressions hybrides
*  - Appeller define_var3 plutot que define_var2
*  - Facteur de compaction entierement traitee dans define_var3
*  - Ajouter la routine auxiliary_apb
*  Bernard Dugas mai 2012 :
*  Creer la variable time_bnds
*  Bernard Dugas fevrier 2009 :
*  Modification a l'appel de affecte_attr
*  Bernard Dugas juillet 2008 :
*  Initialiser cccpack a npack. On utilise toujours cette variable
*  lorsqu'elle est differente de 999 (elle a donc ete lue en entree)
*  Bernard Dugas ete/automne/hiver 2007 :
*  - Ajouter un argument correspondant au numero d'unite I/O PHIS_UNIT
*    et changer le nom de la routine a define_netcdf2
*  - Supporter les reperes geographiques tournes
*  - Supporter les fichiers standards CMC/RPN
*  Guy Bergeron juillet 2004 : Hybrid height coordiante
*  Anne Frigon juillet 2004 : correction de global_attributs pour global_attributes
*
******

      integer nt,ntim,nk,nlev
      integer ilen,nbr,nattrs
      integer i,ii,oldmod,status,nlen
      integer id,varid,dimid,dimids(maxdim),dimidb(2)
      character  cccname*8,mot_cle*60,ip1out*15,timout*14

      character string*512
      integer   idummy
      integer*1 i1dummy
      integer*2 i2dummy
      real*4    rdum
      real*8    ddummy

      character*4 gethic
      external    gethic
      
*-----------------------------------------------------------------------
*     Evaluer les dimensions et les coordonnees :

      call def_dim_coord(funit,ibuf)

*     Creer un fichier netCDF en mode "define"

      status = nf_create(netcdf_file,nf_clobber,ncid)
      call handle_err2(status,'define_netcdf2')      

*     Definir les dimensions dans le fichier netCDF :

      do id=1,ndims
         if (dim(id)%name .eq. 'time') then
            ilen=nf_unlimited
         else
            ilen=dim(id)%len
         endif
         status = nf_def_dim(ncid, dim(id)%name, ilen, dimid)
         call handle_err2(status,'define_netcdf2')
      enddo


*     Definir les coordonnees et leurs attributs dans le fichier netCDF:

      do id=1,ncoord

         status=nf_def_var(ncid, coord(id)%name, coord(id)%type, 
     .                           coord(id)%ndim, coord(id)%dimid, varid)
         call handle_err2(status,'define_netcdf2')

         call attribut_coord(id,coord(id)%nattr)
         
         call put_attribut(ncid,varid,coord(id)%nattr)

      enddo


*     Definir les variables necessaires a "formula terms" et leurs 
*     attributs dans le fichier netCDF.

      nbr=0
      if (level_desc .eq. "Gal-Chen Levels") then

         nbr=nbr+1
         id=nbr                                            ! ID de "eta"

         do i=1,ndims
            if(dim(i)%name.eq.'zt')ii=i
         enddo

         dimids(1)=xdid
         dimids(2)=ydid

         call affecte_var(id,'eta',nf_float,1,zdid,0)

         call affecte_var(id+1,'tau',nf_float,1,zdid,0)

         call affecte_var(id+2,'htoit',nf_float,1,ii,0)

         call affecte_var(id+3,'h0',nf_float,2,dimids,0)


         call hybrid_vertical_coordinate
     . ( phis_unit, ibuf, id,id+1,id+2,id+3 )

         call attribut_var(id,var(id)%nattr,'auxiliary')   ! attributs de "eta"

         status=nf_def_var(ncid, var(id)%name, var(id)%type, 
     .                           var(id)%ndim, var(id)%dimid, varid)
         call handle_err2(status,'define_netcdf2')

         call put_attribut(ncid,varid,var(id)%nattr)
*
         nbr=nbr+1
         id=nbr                                            ! ID de "tau"

         call attribut_var(id,var(id)%nattr,'auxiliary')   ! attributs de "tau"

         status=nf_def_var(ncid, var(id)%name, var(id)%type, 
     .                           var(id)%ndim, var(id)%dimid, varid)
         call handle_err2(status,'define_netcdf2')

         call put_attribut(ncid,varid,var(id)%nattr)
          
         nbr=nbr+1
         id=nbr                                            ! ID de "htoit"

         call attribut_var(id,var(id)%nattr,'auxiliary')   ! attributs de "tau"

         status=nf_def_var(ncid, var(id)%name, var(id)%type, 
     .                           var(id)%ndim, var(id)%dimid, varid)
         call handle_err2(status,'define_netcdf2')

         call put_attribut(ncid,varid,var(id)%nattr)
 
         nbr=nbr+1
         id=nbr                                            ! ID de "h0"

         call attribut_var(id,var(id)%nattr,'auxiliary')   ! attributs de "h0"

         status=nf_def_var(ncid, var(id)%name, var(id)%type, 
     .                           var(id)%ndim, var(id)%dimid, varid)
         call handle_err2(status,'define_netcdf2')

         call put_attribut(ncid,varid,var(id)%nattr)

      else if (level_desc == "Hybrid Levels"
     .    .or. level_desc == "Log Pressure Hybrid Levels") then

         nbr=nbr+1
         id=nbr                                            ! ID de "ap"

         call affecte_var(id,'ap',nf_double,1,zdid,0)

         call attribut_var(id,var(id)%nattr,'char')        ! attributs de "ap"

         status=nf_def_var(ncid, var(id)%name, var(id)%type, 
     .                           var(id)%ndim, var(id)%dimid, varid)
         call handle_err2(status,'define_netcdf2')

         call put_attribut(ncid,varid,var(id)%nattr)

         nbr=nbr+1
         id=nbr                                            ! ID de "b"

         call affecte_var(id,'b',nf_double,1,zdid,0)

         call attribut_var(id,var(id)%nattr,'char')        ! attributs de "b"

         status=nf_def_var(ncid, var(id)%name, var(id)%type, 
     .                           var(id)%ndim, var(id)%dimid, varid)
         call handle_err2(status,'define_netcdf2')

         call put_attribut(ncid,varid,var(id)%nattr)

         call auxiliary_apb( 'ap','b',nbr )

      endif

*     Definir les coordonnees auxiliaires et leurs attributs dans le
*     fichier netCDF :

      if (project%name.eq.'rotated_latitude_longitude' .or.
     .    project%name.eq.'rotated_pole'               .or.
     .    project%name.eq.'polar_stereographic'      ) then

         dimids(1)=xdid
         dimids(2)=ydid

         nbr=nbr+1
         id=nbr                                            ! ID de "lon"
         
         call affecte_var(id,trim(lon),nf_double,2,dimids,0)
         call affecte_var(id+1,trim(lat),nf_double,2,dimids,0)

         call auxiliary_coordinate(dim(xdid)%len,dim(ydid)%len,nbr+1)
         
         call attribut_var(id,var(id)%nattr,'auxiliary')   ! attributs de "lon"

         status=nf_def_var(ncid, var(id)%name, var(id)%type, 
     .                           var(id)%ndim, var(id)%dimid, varid)
         call handle_err2(status,'define_netcdf2')

         call put_attribut(ncid,varid,var(id)%nattr)
*
         nbr=nbr+1
         id=nbr                                            ! ID de "lat"

         call attribut_var(id,var(id)%nattr,'auxiliary')   ! attributs de "lat"

         status=nf_def_var(ncid, var(id)%name, var(id)%type, 
     .                           var(id)%ndim, var(id)%dimid, varid)
         call handle_err2(status,'define_netcdf2')

         call put_attribut(ncid,varid,var(id)%nattr)
*
         nbr=nbr+1
         id=nbr                                     

         call affecte_var(id,project%name,nf_char,0,dimids,0)

         call attribut_var(id,var(id)%nattr,'projection') 

         status=nf_def_var(ncid, var(id)%name, var(id)%type, 
     .                           var(id)%ndim, var(id)%dimid, varid)
         call handle_err2(status,'define_netcdf2')

         call put_attribut(ncid,varid,var(id)%nattr)

      endif


*     Define the time_bnds variable :

      if (time_bnds_L) then
            
         dimidb=(/ bndsdid , timedid /)
         status=nf_def_var(ncid, 'time_bnds', nf_double, 
     .                           2, dimidb, varid)
         call handle_err2(status,'define_netcdf2')
         
      endif

*     Definir les variables et leurs attributs dans le fichier netCDF :

      do i=1,maxvar

         if(infvar(i)%var_ok) then

            nbr=nbr+1
            id=nbr
                                 
            range(1) = infvar(i)%range(1)
            range(2) = infvar(i)%range(2)

            nk       = infvar(i)%zid
            nlev     = infvar(i)%len(nk)

            if(nlev        > 1            .or.
     .         infvar(i)%unique_L         .or.
     .         level_desc == "Surface"    .or.
     .         level_desc == "Sea Level"  .or.
     .         ccc_pktyp(1:2) /= 'SQ'   ) then
               mot_cle='dummy'
            else ! Transmettre le niveau unique approprie
               call get_ip1_string( infvar(i)%niv(1),ip1out )
               nlen=len( ip1out )
               do ii=1,nlen
                  if (ip1out(ii:ii) /= ' ') exit
               enddo
               ip1out = ip1out(ii:nlen)
               mot_cle='single_level '//trim( ip1out )
            endif

            nt       = infvar(i)%tid
            ntim     = infvar(i)%len(nt)

            if(ntim == 1) then
               write(timout,'(I14.14)') infvar(i)%pas(1)/100
               if (mot_cle == 'dummy') then
                  mot_cle='single_time '//timout
               else
                  mot_cle=trim( mot_cle )//'; single_time '//timout
               endif
            endif

            call define_var3(id,i) !temporairement, var(id)%name=infvar(i)%name 

            call attribut_var(id,var(id)%nattr,mot_cle)    ! sortie : var(id)%name 
                                                           ! correctement definie

            status=nf_def_var(ncid, var(id)%name, var(id)%type, 
     .                               var(id)%ndim, var(id)%dimid, varid)
            call handle_err2(status,'define_netcdf2')
         
            call put_attribut(ncid,varid,var(id)%nattr)
               
         endif

      enddo

      nvars=nbr


*     Definition des Global Attributes
      
*     NOTA : Le type derive "var(maxvar)" est utilise temporairement,
*            le temps d'ecrire les global attributes.

      call affecte_var(maxvar,'global_attributes',nf_char,0,dimids,0)

      call attribut_var(maxvar,var(maxvar)%nattr,'char') 

      if (ccc_pktyp(1:2).eq.'SQ') then

         do i=1,var(maxvar)%nattr
            if (attr(i)%name.eq.'title') then

               if (meta_title /= ' ') then

                  string = meta_title

               else

                  ! Mettre l'etiquette et le nom du fichier dans 
                  ! l'attribut global 'title'.

                  call getnam( funit,string )
                  string(14:128) = string(1:115)
                  string(13:13)  = ' '

                  string(1:04)   = gethic('ETIK1', ibuf )
                  string(5:08)   = gethic('ETIK2', ibuf )
                  string(9:12)   = gethic('ETIK3', ibuf )

               endif

               nlen = len_trim( string )
               call affecte_attr(i-1,nf_char,'title',nlen,string,
     .                       i1dummy,i2dummy,idummy,rdum,ddummy)
            endif
         enddo

         if (meta_title /= ' ' .and. i > var(maxvar)%nattr) then

            ! L'attribut 'title n'a pas ete trouve dans le fichier
            ! attribut_netcdf.dat. On l'ajoute explicitement et
            ! la valeur non-nulle de meta_title y est deposee

            string = meta_title ; nlen = len_trim( string )
            call affecte_attr(var(maxvar)%nattr,nf_char,'title',nlen,
     .                     string,i1dummy,i2dummy,idummy,rdum,ddummy)

         endif

      endif

      call put_attribut(ncid,nf_global,var(maxvar)%nattr)

*     (Peut-etre) Enlever le mode "fill mode"

*     NOTA: Par defaut, le fichier est initialisé par une valeur (_FillValue)
*           ou son defaut (nf_fill_* c.f. netcdf.inc) fonction du type de la
*           varibale. Dans notre cas, lorsque les valeurs manquantes ne sont
*           pas definie (avec les arguments "-fill_ccc" ou "-mvalue"), nous
*           ecrivons des tableau en entier et le mode "fill mode"
*           est alors redondant.

      if (.not.fill_ccc_def) then
         status = nf_set_fill(ncid,nf_nofill,oldmod)
         call handle_err2(status,'define_netcdf2')
      endif

*     Fermeture du mode "define"

         status = nf_enddef(ncid)
         call handle_err2(status,'define_netcdf2')

*-----------------------------------------------------------------------
 9904 format(a4)
*-----------------------------------------------------------------------
      end
      subroutine auxiliary_apb( noma,nomb,nbr )

      implicit none

      integer       nbr
      character*(*) noma,nomb

      include 'cdf2ccc.h'
      include 'infomem.h'
      include 'varmem.h'

      integer  i,iap,ib

      iap = -1 ; ib = -1

      do i=1,nbr
         if (iap > 0 .and. ib > 0) exit
         if(var(i)%name == noma) then
            iap = i
         elseif (var(i)%name == nomb) then
            ib =  i
         endif
      enddo

      if (iap > 0 .and. ib > 0) then
         do i=1,maxlev ! Inverser l'ordre
            variable(i,iap) = cap(maxlev-i+1)
            variable(i,ib ) = cb (maxlev-i+1)
         enddo
      end if

      return
      end

      
