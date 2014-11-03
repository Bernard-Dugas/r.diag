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
      subroutine def_attr_polar_stereographic(nbr)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'

      integer nbr

******
*
*AUTEUR Guy Bergeron             juin  2003
*
*     Definir les attributs de la variable "polar_stereographic"  
*     (i.e definir les parametres de la projection)    
*
*REVISIONS
*
*  B. Dugas fevrier 2009 :
*    Modification a l'appel de affecte_attr
*  B. Dugas juillet 2008 :
*    Corrections au definitions des attributs (selon A. Frigon)
*    - straight_vertical_longitude_from_pole
*    - false_easting et false_northing
*  B. Dugas mai 2007 :
*    Ajout de latproj (latitude_of_projection_origin)
*  A. Frigon Juin 2006 : 
*    Ajout de descripteurs de grille polaire stereo : NHEM et D60
*    sous les noms hemisphere_of_standard_parallel et resolution_at_standard_parallel
*    car vers netcdf tout est general pour PS nord/sud 
*    selon IBUF(7) lu et assigne a nhem
*    dans def_dim_coord.f
*
******

      integer*1    i1dummy
      integer*2    i2dummy
      integer      idummy
      real*4       rdum
      real*8       ddummy
      character*80 dummy

      integer      i
      real         lamda0

      real      pi ! distance du pole selon x en nbre de dx
      real      pj ! distance du pole selon y en nbre de dy
      real     d60 ! valeur de dx vrais a 60 degres de l'hemisphere nhem
      real    dgrw ! angle entre l'axe des x et Greenwich(degres ouest positif)
      real latproj ! latitude_of_projection_origin (= either +90. or -90.)
      real    nhem ! hemisphere nord (1), sud (2) ou global (0)
*
*-----------------------------------------------------------------------

      latproj = -1024.

      do i=1,project%len
         if (project%nampar(i).eq.'pi'     ) pi      = project%value(i)
         if (project%nampar(i).eq.'pj'     ) pj      = project%value(i)
         if (project%nampar(i).eq.'d60'    ) d60     = project%value(i)
         if (project%nampar(i).eq.'dgrw'   ) dgrw    = project%value(i)
         if (project%nampar(i).eq.'latproj') latproj = project%value(i)
         if (project%nampar(i).eq.'nhem'   ) nhem    = project%value(i)
      enddo

      lamda0 = -(dgrw+90)
      if (lamda0 < -180) lamda0=lamda0+360
      if (lamda0 >  180) lamda0=lamda0-360

      if (latproj < -1000.) then
         if (nhem == 1) latproj = +90.
         if (nhem == 2) latproj = -90.
      endif

      call affecte_attr(nbr,nf_char,
     .                 'grid_mapping_name',len_trim( project%name ),
     .                  project%name,i1dummy,i2dummy,idummy,rdum,ddummy)

      call affecte_attr(nbr,nf_float,
     .                 'straight_vertical_longitude_from_pole',
     .                  1,dummy,i1dummy,i2dummy,idummy,lamda0,ddummy)

      call affecte_attr(nbr,nf_float,
     .                 'standard_parallel',
     .                  1,dummy,i1dummy,i2dummy,idummy,60.0,ddummy)

      call affecte_attr(nbr,nf_float,
     .                 'false_easting',1,dummy,
     .                  i1dummy,i2dummy,idummy,(pi-1.)*d60,ddummy)

      call affecte_attr(nbr,nf_float,
     .                 'false_northing',1,dummy,
     .                  i1dummy,i2dummy,idummy,(pj-1.)*d60,ddummy)


***debug        write(6,*) "project%value======= pour nhem",nhem           !debug

      call affecte_attr(nbr,nf_float,
     .                 'latitude_of_projection_origin',
     .                  1,dummy,i1dummy,i2dummy,idummy,latproj,ddummy)

      call affecte_attr(nbr,nf_float,
     .                 'resolution_at_standard_parallel',
     .                  1,dummy,i1dummy,i2dummy,idummy,d60,ddummy)

      return
*-----------------------------------------------------------------------

      end
      subroutine def_attr_rotated_lat_lon(nbr)

      implicit none

      include 'netcdf.inc'
      include 'cdf2ccc.h'
      include 'ztypmem.h'

      integer nbr
*
*     Bernard Dugas            mai   2007
*
*     Definir les attributs de la variable "rotated_pole"  
*     (i.e definir les parametres de la projection)
*     
*
*REVISIONS
*

      integer*1    i1dummy
      integer*2    i2dummy
      integer      idummy
      real*4       rdum
      real*8       ddummy
      character*80 dummy

*-----------------------------------------------------------------------
      call affecte_attr(nbr,nf_char,
     .                 'grid_mapping_name',26,
     .                 'rotated_latitude_longitude',
     .                  i1dummy,i2dummy,idummy,rdum,ddummy)

      call affecte_attr(nbr,nf_float,
     .                 'grid_north_pole_latitude',
     .                  1,dummy,i1dummy,i2dummy,idummy,gnplat,ddummy)

      call affecte_attr(nbr,nf_float,
     .                 'grid_north_pole_longitude',
     .                  1,dummy,i1dummy,i2dummy,idummy,gnplon,ddummy)

      call affecte_attr(nbr,nf_float,
     .                 'north_pole_grid_longitude',
     .                  1,dummy,i1dummy,i2dummy,idummy,longpol,ddummy)

      return
*-----------------------------------------------------------------------

      end
