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
      subroutine auxiliary_coordinate(NI,NJ,nbr)

      implicit none

      include 'cdf2ccc.h'
      include 'infomem.h' 
      include 'varmem.h'
      include 'ztypmem.h'

      integer ni,nj,nbr


******
*
*AUTEUR Guy Bergeron         juillet 2003
*
*     Evaluer les variables de coordonnees auxiliaires lon(x,y) et lat(x,y) 
*     associees a xc et yc dans le cas polaire-stereographique ou les tictacs
*     >> et ^^ associes a rlon et rlat dans le cas d'une grille a pole tourne
*
* REVISIONS
*
*     B. Dugas  juillet 2012 :
*     - Lors du calcul de longpol, tenir compte de la rotation de -rlonoff
*       qui est appliquee a la coordonnee xid dans eval_rlonlat
*     - Corriger la boucle ou on force "-180. < lonr < 180."
*     B. Dugas  mai 2009 :
*       Tenir compte de invj dans les definitions de variable(:,lonid)
*       et de variable(:,latid) dans le cas "polar_stereographic"
*     B. Dugas  novembre 2008 :
*       Tenir compte de invj dans les definitions de variable(:,lonid)
*       et de variable(:,latid) dans le cas "rotated_pole"
*     B. Dugas  automne 2007 :
*     - Faire un "call xit" en cas d'erreur
*     - Ajouter le support des grille a repaires geographiques tournes
*     - Tenir compte de la variable latproj pour champs polaire-stereographiques
*     A. Frigon juin 2006 : 
*       Corrige commentaire d60 "vrais a 60N" par "vrais a 60 deg de l'hemisphere nhem"
*       car vers netcdf tout est general pour PS nord/sud 
*       selon IBUF(7) lu et assigne a nhem
*       dans def_dim_coord.f
*
******

      integer i,j,ij,iji,jinv,lonid,latid
      real    is,js,xdx,ydy,rdummy,xlat,xlon

      real      pi ! distance du pole selon x en nbre de dx
      real      pj ! distance du pole selon y en nbre de dy
      real     d60 ! valeur de dx vrais a 60 degres de l'hemisphere nhem
      real    dgrw ! angle entre l'axe des x et Greemwich(degres ouest positif)
      real latproj ! latitude_of_projection_origin (= either +90. or -90.)
      integer nhem ! hemisphere nord (1), sud (2) ou global (0)
      integer  nis ! nbre de points en X grille de type f
      integer  njs ! nbre de points en Y grille de type f

*     coordonnees vraies de deux points sur l'equateur tourne.
*     Le premier est au centre de la grille tournee et le second
*     ailleurs sur cet equateur
      real dlon1,dlat1 ! true longitude,latitude of first point
      real dlon2,dlat2 ! true longitude,latitude of second point
      real*8 xyz(3)
      real   latgpol

      integer LIJ
      LIJ(i,j)=(j-1)*ni+i

      LOGICAL              rpn_info
      COMMON     /ZZVERBO/ rpn_info
*-----------------------------------------------------------------------
      lonid=-999

      do i=1,project%len

         if (project%nampar(i).eq.'nis'    ) nis     = project%value(i)
         if (project%nampar(i).eq.'njs'    ) njs     = project%value(i)
         if (project%nampar(i).eq.'pi'     ) pi      = project%value(i)
         if (project%nampar(i).eq.'pj'     ) pj      = project%value(i)
         if (project%nampar(i).eq.'d60'    ) d60     = project%value(i)
         if (project%nampar(i).eq.'dgrw'   ) dgrw    = project%value(i)
         if (project%nampar(i).eq.'latproj') latproj = project%value(i)
         if (project%nampar(i).eq.'nhem'   ) nhem    = project%value(i)

         if (project%nampar(i).eq.'dlon1'  ) dlon1   = project%value(i)
         if (project%nampar(i).eq.'dlat1'  ) dlat1   = project%value(i)
         if (project%nampar(i).eq.'dlon2'  ) dlon2   = project%value(i)
         if (project%nampar(i).eq.'dlat2'  ) dlat2   = project%value(i)

      enddo
*
      do i=1,nbr
         if(var(i)%name.eq.'lon' .or.
     .      var(i)%name.eq.'LON') then
            lonid=i
         else if(var(i)%name.eq.'lat' .or.
     .           var(i)%name.eq.'LAT') then
            latid=i 
         endif
      enddo
            
      if (lonid.eq.-999) then
         print *,' Probleme avec auxilary_coordinate '
         call xit('auxiliary_coordinate',-1)
      endif
*
      if (project%name.eq.'rotated_latitude_longitude'  .or.
     .    project%name.eq.'rotated_pole'              ) then

         if (RPN_INFO)
     .   write(6,6000) dlon1,dlat1,dlon2,dlat2

*        au retour de d_rota, rrot contient la matrice
*        de rotation cartesienne et gnplon,gnpla sont
*        definis

         call d_rota( lonr,latr, alon,alat, 
     +                dlon1,dlat1,dlon2,dlat2,
     +                gnplon,gnplat, ni,nj )

         if (abs( gnplat ) < 1.e-6 ) gnplat = 0.0

         do ij=1,ni*nj
            if (lonr(ij)-180. > 0.001) lonr(ij)=lonr(ij)-360.0
            if (lonr(ij)+180. <-0.001) lonr(ij)=lonr(ij)+360.0
         enddo

         ! Tenir compte du "F_lon(i) = amod( F_lon(i) , 360.0 )"
         ! qui est effectue dans d_cartall (appelle par d_rota)

         do j=1,nj
            if (lonr(LIJ(ni,j))   < lonr(LIJ(ni-1,j)) .and.
     +     abs( lonr(LIJ(ni,j)) ) < 0.001) 
     +          lonr(LIJ(ni,j))   = 360.  
         enddo
         

*        xyz contiendra la troisieme collonne de la matrice
*        de rotation rrot et puisque

*        rrot: (repere non-tourne) --> (repere tourne)

*        alors 
*        rrot x T(0,0,1) = coordonnes cartesiennes du pole
*                          non-tourne dans le repere tourne

         xyz(:) = rrot(:,3)
         call d_cartall( longpol,latgpol,xyz,1 )

         ! Tenir compte de la rotation de -rlonoff
         ! qui est appliquee a la coordonnee dans
         ! la routine eval_rlonlat

         longpol = longpol+rlonoff
         if (longpol-360.+rlonoff > 0.001) longpol=longpol-360.0
         if (longpol     +rlonoff <-0.001) longpol=longpol+360.0

         ij=0
         do j=1,nj
            if (invj) then
               iji=(nj-j)*ni    ! inversion des y
            else
               iji=(j-1)*ni
            endif
            do i=1,ni
               ij=ij+1
               iji=iji+1
               variable(iji,lonid)=lonr(ij)
               variable(iji,latid)=latr(ij)
            enddo
         enddo

      elseif (project%name.eq.'polar_stereographic') then

         is=(nis-ni)*0.5
         js=(njs-nj)*0.5

         if (latproj.eq.90.) then
            nhem = 1
         else if (latproj.eq.-90.) then
            nhem = 2
         endif
      
         do i=1,ni
         do j=1,nj
            if (invj) then            
               ij=i+(nj-j)*ni   ! inversion des y
            else
               ij=i+(j-1)*ni
            endif

            xdx=real(i)+is-pi
            ydy=real(j)+js-pj            
            call llfxy(xlat,xlon, xdx,ydy,d60,dgrw,nhem)
     .                                          
            variable(ij,latid) = xlat
            variable(ij,lonid) = xlon

            if (variable(ij,lonid).lt.0.) 
     .                        variable(ij,lonid)=variable(ij,lonid)+360.
            
         enddo
         enddo

      else
         write(6,6001) project%name
      endif

*-----------------------------------------------------------------------
 6000 format(/21x,'dlon1 = ',f10.2/
     .        21x,'dlat1 = ',f10.2/
     .        21x,'dlon2 = ',f10.2/
     .        21x,'dlat2 = ',f10.2/)
 6001 format('auxiliary_coordinate : grille non definie pour grid_desc '
     .                                                             ,a20)
*-----------------------------------------------------------------------
      end
