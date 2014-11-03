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
      subroutine get_dim(NCID)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'dimmem.h'
      
      integer ncid

******
*
*AUTEUR Guy Bergeron    Juin 2003
*    
*     Affecte les valeurs au type derive dim(i) a partir du fichier netCDF.
*
******

      integer status,dimid,len
      character*80 name

*-----------------------------------------------------------------------

      do dimid=1,ndims
         status = nf_inq_dim(NCID,DIMID,name,len) ! get dimname
         call handle_err2(status,'get_dim')
         if(name .eq. 'num_values') then 
            spec=.true.
            numdid=dimid
            invj=.false.           ! pas d'inversion  pour cas spec
         endif
         dim(dimid)%name=name                     ! nom de la dimension 
         dim(dimid)%len=len                       ! longueur de la dimension
      enddo
*-----------------------------------------------------------------------
      end
      subroutine test_dim ()

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'
      
* AUTEUR Bernard Dugas    Juin 2008
*     
*     tenter de definir les dimensions horizontale si 
*     la routine get_coord2 n'y est pas arrive. Les
*     noms des dimensions sont alors utilises
*
* REVISIONS
*
*  B.Dugas aout '12 :
*  - Ajouter 'height' a la liste des coordonnees verticale reconnues
*  B.Dugas juillet '12 :
*  - Tenir compte de unlimdimid dans la recherche des dimensions:
*    Seul le descripteur TID peut etre associe a cette dimension
*  B.Dugas mai '12 :
*  - Reconnaitre 'plev' comme nom possible la coordonne verticale
*  B.Dugas mars '10 :
*  - Reconnaitre 't' comme nom possible la coordonne temporelle.
*  B.Dugas mai '09 :
*  - Coordonnees NetCDF "[xyzt]coord" specifiees en arguments ?
*  B.Dugas oct '08 :
*  - Supporter la coordonnee verticale ayant pour nom "lev"
*
******
      integer dimid,nlen,tdid,ncc
      character(len=128) cfield

*-----------------------------------------------------------------------

      if (coord(xid)%dimid(1).ne.-1 .and.
     .    coord(yid)%dimid(1).ne.-1 .and.
     .    coord(zid)%dimid(1).ne.-1 .and.
     .    coord(tid)%dimid(1).ne.-1) then

!        On a tout trouve

         return

      else

         ncc = ncoord

         do dimid=1,ndims

            call clean_char( dim(dimid)%name,cfield,nlen )
            call up2low( cfield, cfield )

!           Chercher coordonnee T (avec caracteristique unlimdimid)

            if (coord(tid)%dimid(1).eq.-1) then

               if((cfield   == 'time'         .or. 
     .             cfield   == 't'            .or.
     .             cfield   ==  tcoord)       .and.
     .             dimid    ==  unlimdimid)   then

                  ncc = ncc+1
                  tid = ncc

                  coord(tid)%ndim     = 1
                  coord(tid)%nattr    =-1
                  coord(tid)%dimid(1) = dimid

                  timedid             =     dimid
                  coord(tid)%name     = dim(dimid)%name

                  cycle

               endif

            endif

!           Chercher coordonnee X

            if (coord(xid)%dimid(1).eq.-1) then

               if((cfield == 'x'          .or.
     .             cfield == 'x_2'        .or.
     .             cfield == 'lon'        .or.
     .             cfield == 'rlon'       .or.
     .             cfield == 'longitude'  .or. 
     .             cfield == 'west_east'  .or. 
     .             cfield ==  xcoord    ) .and.
     .             dimid  /=  unlimdimid) then

                  ncc = ncc+1
                  xid = ncc

                  coord(xid)%ndim     = 1
                  coord(xid)%type     = nf_float
                  coord(xid)%nattr    =-1
                  coord(xid)%mult     = 1.0
                  coord(xid)%add      = 0.0
                  coord(xid)%dimid(1) = dimid

                  coord(xid)%name     = dim(dimid)%name
                  lon                 = cfield

                  cycle

               endif

            endif

!           Chercher coordonnee Y

            if (coord(yid)%dimid(1).eq.-1) then

               if((cfield == 'y'            .or.
     .             cfield == 'y_2'          .or.
     .             cfield == 'lat'          .or.
     .             cfield == 'rlat'         .or.
     .             cfield == 'latitude'     .or. 
     .             cfield == 'south_north'  .or. 
     .             cfield ==  ycoord      ) .and.
     .             dimid  /=  unlimdimid  ) then

                  ncc = ncc+1
                  yid = ncc

                  coord(yid)%ndim     = 1
                  coord(yid)%type     = nf_float
                  coord(yid)%nattr    =-1
                  coord(yid)%mult     = 1.0
                  coord(yid)%add      = 0.0
                  coord(yid)%dimid(1) = dimid

                  coord(yid)%name     = dim(dimid)%name
                  lat                 = cfield

                  cycle

               endif

            endif

!           Chercher coordonnee Z

            if (coord(zid)%dimid(1).eq.-1) then

               if((cfield == 'p'           .or.
     .             cfield == 'lev'         .or.
     .             cfield == 'plev'        .or.
     .             cfield == 'level'       .or.
     .             cfield == 'height'      .or.
     .             cfield == 'bottom_top'  .or.
     .             cfield ==  zcoord     ) .and.
     .             dimid  /=  unlimdimid ) then

                  ncc = ncc+1
                  zid = ncc

                  coord(zid)%ndim     = 1
                  coord(zid)%type     = nf_float
                  coord(zid)%nattr    =-1            
                  coord(zid)%mult     = 1.0
                  coord(zid)%add      = 0.0
                  coord(zid)%dimid(1) = dimid

                  zdid                =     dimid
                  coord(zid)%name     = dim(dimid)%name

                  cycle

               endif

            endif

         enddo

         if (coord(zid)%dimid(1) == -1) then

            ! initialisation pour cas a 1 niveau
            
            zdid=maxdim

            dim(zdid)%name='lev'
            dim(zdid)%len=1
            dim(zdid)%duplic=0

            coord(zid)%name    = dim(zdid)%name
            coord(zid)%type    = nf_float
            coord(zid)%nattr   = 0            
            coord(zid)%ndim    = 1
            coord(zid)%mult    = 1.0
            coord(zid)%add     = 0.0
            coord(zid)%dimid(1)= zdid

         end if

         if (coord(tid)%dimid(1) == -1) then

            ! initialisation pour cas intemporel

            tdid=maxdim-3
            timedid=tdid

            dim(tdid)%name='time'
            dim(tdid)%len=1
            dim(tdid)%duplic=0

            coord(tid)%name    = dim(tdid)%name
            coord(tid)%type    = nf_float
            coord(tid)%nattr   = 0            
            coord(tid)%ndim    = 1
            coord(tid)%mult    = 1.0
            coord(tid)%add     = 0.0
            coord(tid)%dimid(1)= tdid

         endif

      endif

      return
      end
