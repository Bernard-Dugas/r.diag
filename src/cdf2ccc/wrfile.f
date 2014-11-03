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
      subroutine wrfile(NCID,funit)  !,cccname)
      
      implicit none


      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h' 
      include 'varmem.h'
      include 'ibufmem.h'
      include 'workmem.h'     !lag 

      integer ncid,funit


******
*
*AUTEUR Guy Bergeron         juin 2003
*
*     Ecriture des valeurs des variables dans le fichier netCDF
*
*REVISIONS 
*
*  Bernard Dugas juin 2013 :
*  - Appeller des versions de lire_cccvar, put_vard et put_varda qui
*    tiennent compte de fill_ccc_def dans leurs traitements de variable
*  Bernard Dugas juillet 2012 :
*  - Corriger les definitions de start et count, en tenant
*    compte d'un (possiblement) nombre variable de dimensions
*     - Tenir compte des variables intemporelles "meteorologiques"
*    qui doivent tout de meme etre lues sur le fichier funit.
*    Ces variables incluent maintenant toutes celles qui
*    ne sont definies qu'a un seul moment.
*  Bernard Dugas juin 2012 :
*  - Tenir compte du nombre de niveaux verticaux nlev
*    dans la mesure ou lorsque nlev <= 0, la verticale
*    n'est plus declaree comme dimension
*  - Appeller lire_cccvar2 plutot que lire_cccvar
*  Bernard Dugas mai 2012 :
*  - Sauver time_bnds et re-definir la variables
*    dcoordonne(:,tid) lorsque time_bnds_L est vrai
*  - Deplacer l'ecriture de 'time' et 'time_bnds'
*  Bernard Dugas fevrier 2009 :
*  - CONV est vrai pour les donnees nf_byte
*  Bernard Dugas avril 2008 :
*  - Tenir compte des cas ou les attributs scale_factor
*    et add_offset sont en format "real*8"
*  Bernard Dugas juillet 2007 :
*  - On utlise attunit plutot que "99"
*  - Allocation de memoire de type "automatic array"
*  Guy Bergeron juillet 2004 : Restructuration de l'algorithme d'ecriture     
*
******

******netCDF

      integer i,j,itime,nlen,indx1(2),indx2(2),nlev
      integer status,id,varid,timeid,bndsid,vartid
      integer start(maxdim),count(maxdim)

      real*8  scale,offset,dum1,dum2, fillval
      real*8  dnewtime       !lag

      character*4 cccname
      character*128 dummy

      logical conv,ok(2),mss

*-----------------------------------------------------------------------
      do i=1,nlist
         if(list(i)%name.eq.'time') then
            list(i)%var_ok=.false.
            timeid=i                ! le id de "time" dans le fichhier netcdf
         else if (list(i)%name.eq.'time_bnds') then
            list(i)%var_ok=.false.
            bndsid=i
         endif
      enddo

*     Sauver time et optionellement tim_bnds

      do itime=1,dim(timedid)%len ! boucle temporelle
         
         dnewtime=dcoordonne(itime,tid)

         if (time_bnds_L) then  ! Implique notamment que tlbl=T

*           Re-Definir la variable de sortie temporelle

            dnewtime = 0.5 * ( time_bnds(2,itime) + time_bnds(1,itime) )

*           Ecrire les bornes temporelles correspondantes

            indx1  = (/ 1 , itime /) ; indx2 = (/ 2 , itime /)
            status = nf_put_var1_double( ncid,bndsid,indx1,
     .                                   time_bnds(1,itime) )
            call handle_err2(status,'wrfile')

            status = nf_put_var1_double( ncid,bndsid,indx2,
     .                                   time_bnds(2,itime) )
            call handle_err2(status,'wrfile')

         endif
                  
         if (.not.leap .and.
     .       .not.tlbl      )call enleve_bissextile2(ncid,
     .                               dcoordonne(itime,tid),dnewtime)       

         status=nf_put_var1_double(ncid,timeid,itime,dnewtime)
         call handle_err2(status,'wrfile')

      enddo

*     1a) Variables de Coordonnee (excepte "time") : 

      scale =1.0
      offset=0.0
      conv  =.false.
      mss   =.false.

      do id=1,ncoord
         if(coord(id)%name.ne.'time') then

            do i=1,nlist
               if(coord(id)%name.eq.list(i)%name) varid=i
            enddo
            list(varid)%var_ok =.false.
            nlen=dim(coord(id)%dimid(1))%len
            call put_vard2(ncid,varid,coord(id)%type,nlen,
     .             dcoordonne(1,id),scale,offset,conv,mss,fillval)
         endif
      enddo

*     1b) Variables intemporelles :

      do id=1,nvars
      
         scale =1.0
         offset=0.0
         conv  =.false.

         do i=1,nlist
            if(var(id)%name.eq.list(i)%name) varid=i 
         enddo

         nlen=1
         do i=1,var(id)%ndim
            if(var(id)%dimid(i).ne. timedid)then
               nlen=nlen*dim(var(id)%dimid(i))%len
            else
               goto 100          ! ne pas traiter les variables temporelles
            endif           
         enddo

         list(varid)%var_ok=.false.

         call get_attribut(ncid,varid,var(id)%nattr,var(id)%name)

         do i=1,var(id)%nattr
            if        ( attr(i)%name.eq.'scale_factor' ) then
               if     ( attr(i)%type.eq.nf_float       ) then
                  scale=attr(i)%rvalue(1)
               else if( attr(i)%type.eq.nf_double      ) then
                  scale=attr(i)%dvalue(1)
               endif
            else if   ( attr(i)%name.eq.'add_offset'   ) then
               if     ( attr(i)%type.eq.nf_float       ) then
                 offset=attr(i)%rvalue(1)
               else if( attr(i)%type.eq.nf_double      ) then
                 offset=attr(i)%dvalue(1)
               endif
            else if   ( attr(i)%name.eq.'_FillValue'   ) then
               if     ( attr(i)%type.eq.nf_double      ) then
                 fillval=attr(i)%dvalue(1)
               else if( attr(i)%type.eq.nf_float       ) then
                 fillval=attr(i)%rvalue(1)
               else if( attr(i)%type.eq.nf_short       ) then
                 fillval=attr(i)%i2value(1)
               else if( attr(i)%type.eq.nf_byte        ) then
                 fillval=attr(i)%i1value(1)
               endif
            endif
         enddo
         
         if(var(id)%type.eq.nf_byte   .or.
     .      var(id)%type.eq.nf_short) conv=.true.

         rewind attunit
         call get_name(attunit,var(id)%name,dummy,cccname,
     .                 dum1,dum2,ok)

         if (ok(1)) then

            call leadblk( cccname ) ; mss = .false.

            nlev=0
            do i=1,var(id)%ndim
               if (var(id)%dimid(i) == zdid) then
                  nlev = dim(var(id)%dimid(i))%len
                  exit
               endif
            enddo

            do i=1,maxvar ! Doit-on lire cette variable ?

               if (infvar(i)%name == cccname) then ! Oui, on lit funit

                  call precede( funit, -1 )

                  itime=1
                  call lire_cccvar3( ncid,funit,id,cccname,
     .                               nlev,itime,ibuf,mss )

                  exit

               endif

            enddo
         endif
               
         if(var(id)%ndim > 0)
     .   call put_vard2(ncid,varid,var(id)%type,nlen,variable(1,id),
     .                                scale,offset,conv,mss,fillval)

 100  enddo


*     2-Variables meteorologiques:

      do varid=1,nlist

         if(list(varid)%var_ok) then
            
            scale =1.0
            offset=0.0
            conv  =.false.

            do i=1,nvars
               if(var(i)%name.eq.list(varid)%name) id=i 
            enddo

            write(6,6000)trim(var(id)%name)
            
            list(varid)%var_ok=.false.

            nlen=1
            do i=1,var(id)%ndim
               if(var(id)%dimid(i).ne. timedid)
     .                               nlen=nlen*dim(var(id)%dimid(i))%len
            enddo

            nlev=0
            do i=1,var(id)%ndim
               if (var(id)%dimid(i) == zdid) then
                  nlev = dim(var(id)%dimid(i))%len
                  exit
               endif
            enddo

            call get_attribut(ncid,varid,var(id)%nattr,var(id)%name)
            
            do i=1,var(id)%nattr
               if        ( attr(i)%name.eq.'scale_factor' ) then
                  if     ( attr(i)%type.eq.nf_float       ) then
                     scale=attr(i)%rvalue(1)
                  else if( attr(i)%type.eq.nf_double      ) then
                     scale=attr(i)%dvalue(1)
                  endif
               else if   ( attr(i)%name.eq.'add_offset'   ) then
                  if     ( attr(i)%type.eq.nf_float       ) then
                    offset=attr(i)%rvalue(1)
                  else if( attr(i)%type.eq.nf_double      ) then
                    offset=attr(i)%dvalue(1)
                  endif
               else if   ( attr(i)%name.eq.'_FillValue'   ) then
                  if     ( attr(i)%type.eq.nf_double      ) then
                     fillval=attr(i)%dvalue(1)
                  else if( attr(i)%type.eq.nf_float       ) then
                     fillval=attr(i)%rvalue(1)
                  else if( attr(i)%type.eq.nf_short       ) then
                     fillval=attr(i)%i2value(1)
                  else if( attr(i)%type.eq.nf_byte        ) then
                     fillval=attr(i)%i1value(1)
                  endif
               endif
            enddo
         
            if(var(id)%type.eq.nf_byte   .or.
     .         var(id)%type.eq.nf_short) conv=.true.

*           Initialiser start et count:

            vartid = -1
            do i=1,ndims
               do j=1,var(id)%ndim
                  if (var(id)%dimid(j) == coord(i)%dimid(1)) then
                     if (var(id)%dimid(j) == timedid) vartid=j
                     count(j)=dim(i)%len ! longueurs des boucles
                     start(j)=1 ! indice du debut de l'affectation
                     cycle
                  endif
               enddo
            enddo

            if (vartid > 0) then
               count(vartid)=1
            else
               call xit(' wrfile ',-1 )
            endif

            rewind attunit
            call get_name(attunit,var(id)%name,dummy,cccname,
     .                    dum1,dum2,ok)

            call precede( funit, -1 ) ; mss = .false.

            boucle_temporelle : do itime=1,dim(timedid)%len
         
* NOTA : Il ne faut pas modifier "dcoordonne(i,tid)" p.c.q, utilisee par 
*        recget dans lire_cccvar

               start(vartid)=itime

               call lire_cccvar3( ncid,funit,id,cccname,
     .                            nlev,itime,ibuf,mss)
         
*        Definir la variable i2val = int((value - offset)*scale) et 
*        ecriture dans le fichier de sortie (lorsque npack = -16).

               call put_varda2(ncid,varid,var(id)%type,nlen,start,
     .            count,variable(1,id),scale,offset,conv,mss,fillval)

            enddo boucle_temporelle

         endif
      enddo
*-----------------------------------------------------------------------

 6000 format(/,' VARIABLES : ',a)

*-----------------------------------------------------------------------
      end
