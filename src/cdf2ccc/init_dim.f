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
      subroutine init_dim

      implicit none

      include 'cdf2ccc.h'
      include 'dimmem.h'
      include 'infomem.h'

******
*
* AUTEUR Guy Bergeron    Juin 2003
*      
*     initialisation par defaut de certaines variables.
*
* REVISIONS
*
*  Bernard Dugas Aout 2018 :
*   Ne plus initialiser les variables project%nampar(project%len) et
*   project%value(project%len) que si project%name = 'unknown'     
*
*  Bernard Dugas Octobre 2007 :
*   Initialiser coord(xid)%dimid(1)et coord(yid)%dimid(1) a -1
*
******

      integer i,tdid
*-----------------------------------------------------------------------
*
      do i=1,maxdim
         dim(i)%name='xxx'              ! noms de la dimension 
         dim(i)%len=0                   ! longueur de la dimension
         dim(i)%duplic=0                ! repetition en bout de ligne
      enddo

      tid=maxdim-3
      xid=maxdim-2
      yid=maxdim-1
      zid=maxdim                        

      coord(tid)%dimid(1)= -1
      coord(zid)%dimid(1)= -1
      coord(xid)%dimid(1)= -1
      coord(yid)%dimid(1)= -1

      if (project%name == 'unknown') then
         project%nampar(project%len)='nhem'
         project%value(project%len)=float(0) ! representation globale (defaut)
      endif

*-----------------------------------------------------------------------
      end
