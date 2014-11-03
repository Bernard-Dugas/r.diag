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
      subroutine eval_lonlat(ilon,ilat)

      implicit none

      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'
      include 'varmem.h'

      integer ilat,ilon

******
*
*AUTEUR Guy Bergeron   juin 2003
*
*     Evalue les longitudes et les latitudes definissant la grille "lon/lat"
*     ou "gaussian"
*
*REVISIONS
*
*     Bernard Dugas juillet 2012 :
*     - Limiter la modif precedente aux grilles tournees puisque
*       autrement cela produit une variable coordonnee non monotone
*     Bernard Dugas mai 2012 :
*     - Forcer les longitudes a etre entre -180 et +180
*     Bernard Dugas autome 2007 :
*     - Allocate/deallocate plutot que hpalloc/hpdeallc
*     - On distingue entre les grilles CMC/RPN de type 'A' ou 'B'
*     Bernard Dugas 3 aout 2007 : Tenir compte de la cle invj
*     G. Bergeron ???? : Grille lon/lat regionale
*
******

      real*8 shftlg,shftlt
      data shftlg /0.0/

      integer i,ilath,nlen,nblen
      real*8 olon,olat,dlon,dlat

      real*8  val,pi
      real*8, dimension (:), allocatable :: sl,cl,wl,wossl,rad
*-----------------------------------------------------------------------

      pi = 2.0*asin(1D0)

      nlen=dim(coord(yid)%dimid(1))%len

      if (project%name.eq.'gaussian') then

         allocate( sl(nlen),cl(nlen),wl(nlen),wossl(nlen),rad(nlen) )

         ilath=ilat/2
         call gaussg(ilath,sl,wl,cl,rad,wossl)
         call trigl (ilath,sl,wl,cl,rad,wossl)

         if (invj) then
            do i=1,ilat
               dcoordonne(ilat-i+1,yid)=rad(i)*180./pi ! degrees_north
            end do
         else
            do i=1,ilat
               dcoordonne(     i  ,yid)=rad(i)*180./pi ! degrees_north
            end do
         endif

         do i=1,ilon                               
            dcoordonne(i,xid)=(i-1)*(360.d0/ilon)           ! degrees_east
         end do

         deallocate( sl,cl,wl,wossl,rad )

      else if (project%name(1:14) .eq. 'lon/lat global') then     

         dlon=360d0/ilon

         if (project%name(16:16).eq.'A') then
            dlat=180d0/ilat
            shftlt=0.5*dlat
         else
            dlat=180./float(ilat-1) 
            shftlt=0.0
         endif

         if (invj) then
            do i=1,ilat
               dcoordonne(ilat-i+1,yid)=dble(i-1)*dlat + shftlt-90.0 ! degrees_north
            end do
         else
            do i=1,ilat
               dcoordonne(     i  ,yid)=dble(i-1)*dlat + shftlt-90.0 ! degrees_north
            end do
         endif

         do i=1,ilon
            dcoordonne(i,xid)=dble(i-1)*dlon + shftlg              ! degrees_east
         end do

      else if (project%name .eq. 'lon/lat regional') then         


         do i=1,project%len
            if(project%nampar(i).eq."0lon")olon=project%value(i)
            if(project%nampar(i).eq."0lat")olat=project%value(i)
            if(project%nampar(i).eq."dlon")dlon=project%value(i)
            if(project%nampar(i).eq."dlat")dlat=project%value(i)
         enddo

         do i=1,ilon
            val=dble(i-1)*dlon + olon
            dcoordonne(i,xid)=val                    ! degrees_east
         end do
         
         if (val.gt.360.) then
            do i=1,ilon
               dcoordonne(i,xid)=dcoordonne(i,xid)-360.
            end do
         endif


         do i=1,ilat

            val=dble(i-1)*dlat + olat
            if(val.gt. 90.0)val= 90.0 - (val-90.0)    ! Des valeurs comprises 
            if(val.lt.-90.0)val=-90.0 - (val+90.0)    ! [-90.0,90.0].
            if (invj) then                            ! degrees_north
               dcoordonne(ilat-i+1,yid)=val   
            else
               dcoordonne(     i  ,yid)=val   
            endif

         end do
*     
      else
         write(6,6001) trim(project%name)
         call                                      xit('eval_lonlat',-1)
      end if

!     do i=1,ilon
!        val = dcoordonne(i,xid)
!        if(val > 180.0)val=val-360.0
!        if(val <-180.0)val=val+360.0
!        dcoordonne(i,xid) = val
!     enddo

*-----------------------------------------------------------------------
 6001 format('eval_lonlat : probleme avec grid_desc =-',a,'-')
*-----------------------------------------------------------------------
      end
      subroutine eval_rlonlat(NI,NJ)

      implicit none

      include 'cdf2ccc.h'
      include 'infomem.h'
      include 'varmem.h'
      include 'ztypmem.h'

      integer ni,nj


*     Bernard Dugas  mai 2007
*
*     Evalue les valeurs des coordonnees x(ni) et y(nj) de la projection 
*     definie par project%name (grid_desc) et les descripteurs de grille.
*

* REVISIONS
*     Bernard Dugas juillet 2012 :
*     - Forcer la coordonnee xid a etre entre -rlonoff et 360-rlonoff
*     - La valeur par defaut de rlonoff=-1000. est traduite a 180.
*     Bernard Dugas mai 2012 :
*     - Forcer les longitudes a etre entre -180 et +180

*
******
      integer i,j
      real(8) val
*
*-----------------------------------------------------------------------

      ! La valeur par defaut de rlonoff pour des
      ! grilles tournees est de 180. degres

      if (rlonoff < -999.999) rlonoff = 180.

      ! Appliquer une rotation de -rlonoff en longitudes. La
      ! grille sera alors centree a 180.-rlonoff et les limites
      ! des longitudes seront alors -rlonoff et 360-rlonoff

      do i=1,ni
         dcoordonne(i,xid)=alon(i)-rlonoff
      enddo

      do i=1,ni
         val = dcoordonne(i,xid)
         if (val-360.+rlonoff > 0.001_8) val = val-360.0
         if (val     +rlonoff <-0.001_8) val = val+360.0
         dcoordonne(i,xid) = val
      enddo

      if (invj) then
         do j=1,nj
            dcoordonne(j,yid)=alat(nj-j+1)
         enddo
      else
         do j=1,nj
            dcoordonne(j,yid)=alat(j)
         enddo
      endif

*-----------------------------------------------------------------------
      end
