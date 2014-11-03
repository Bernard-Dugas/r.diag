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
      subroutine eval_xcyc(NI,NJ)

      implicit none

      include 'cdf2ccc.h'
      include 'infomem.h'
      include 'varmem.h'

      integer ni,nj

******
*
* AUTEUR Guy Bergeron   juin 2003
*
*     Evalue les valeurs des coordonnees x(ni) et y(nj) de la projection 
*     definie par project%name (grid_desc) et les descripteurs de grille.
*
* REVISIONS
*
* A. Frigon juin 2006 : 
*   Corrige commentaire d60 "vrais a 60N" par "vrais a 60 deg de l'hemisphere nhem"
*   car vers netcdf tout est general pour PS nord/sud 
*   selon IBUF(7) lu et assigne a nhem
*   dans def_dim_coord.f
*
******

      integer i
      real    is,js       !staggered

      real    d60  ! valeur de dx vrais a 60 deg de l'hemisphere nhem
      integer nis  ! nbre de points en X grille de type f
      integer njs  ! nbre de points en Y grille de type f
*
*-----------------------------------------------------------------------
      do i=1,project%len
         if(project%nampar(i).eq.'nis' )nis =project%value(i)
         if(project%nampar(i).eq.'njs' )njs =project%value(i)
         if(project%nampar(i).eq.'d60' )d60 =project%value(i)
      enddo
*
      is=(nis-ni)*0.5
      js=(njs-nj)*0.5

      do  i=1,ni
         dcoordonne(i,xid)=(i-1+is)*d60
      end do
*     
      do i=1,nj
         dcoordonne(nj-i+1,yid)=(i-1+js)*d60    !inversion des y!!!
c         dcoordonne(i,yid)=(i-1+js)*d60
      end do

*-----------------------------------------------------------------------
      end
