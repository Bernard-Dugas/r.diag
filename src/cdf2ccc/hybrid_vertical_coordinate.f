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
      subroutine hybrid_vertical_coordinate
     .         (phis_unit,ibuf,eta,tau,ht,h0)

      implicit none

      include 'cdf2ccc.h'
      include 'infomem.h' 
      include 'dimmem.h'
      include 'varmem.h'

      integer phis_unit,ibuf(1)

      integer eta,     ! id definissant la variable eta
     .        tau,     ! id definissant la variable tau
     .        ht ,     ! id definissant la variable htoit
     .        h0       ! id definissant la variable h0

******
*
* AUTEUR Guy Bergeron         juin 2004
*
*     Defini les termes necessaires pour "formula terms" associe a 
*     atmospheric hybride height coordinate.
*
*     Sachant que la hauteur geometrique (z) associe a une hauteur 
*     Gal-Chen (ZGC) est definie en fonction de la topgraphie (h0)
*     et de la hauteur du toit (htoit) de la facon suivante :
*
*           z(i,j,k) = h0(i,j) + (1-h0(i,j)/htoit)*ZGC(k)
*
*     il nous est possible de reecrire comme suit :
*     
*           z(i,j,k) = (1-ZGC(k)/htoit)*h0(i,j) + htoit*(ZGC(k)/htoit)
*
*     ou encore
*
*           z(i,j,k) = tau(k)*h0(i,j) + htoit*eta(k)
*
*     ou "tau" et "eta" sont deux termes addimensionnels.
*
*     Pour plus d'information, se referer a :
*     NetCDF Climate and Forecast (CF) Metadata Conventions, appendice D
*     www.cgd.ucar.edu/cms/eaton/cf-metadata/CF-1.0.html 
*
*
* REVISIONS
*
*  Bernard Dugas fevrier 2008 :
*  La lecture du champs de montagne se fait maintenant
*  sur l'unite I/O phis_unit plutot que funit
*
******

      integer k
*

*-----------------------------------------------------------------------
*     evaluer la variable eta et tau:

      do k=1,dim(coord(zid)%dimid(1))%len

         variable(k,eta) = dcoordonne(k,zid)/htoit
         variable(k,tau) = 1.0-variable(k,eta)
      enddo
      
      variable(1,ht) = htoit

*     lire le champs de montagne sur phis_unit

      call get_topo( phis_unit,ibuf, h0 )

*-----------------------------------------------------------------------
      end
